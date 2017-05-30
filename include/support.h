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
#include <cmath>

/** \defgroup tech Technical stuff
*\brief Types, constants, data structs etc
*@{*/

/** \defgroup type Type selection
*@{ 
* \brief Handles type selection to match data
*/
//-----------TYPE HANDLING----------------------------
#ifdef USE_FLOAT

#define ADD_FFTW(x) fftwf_ ## x /**< Add the correct FFTW library prefix to function/variable names */
#define cplx_type ADD_FFTW(complex) /**< Suitable complex type for FFTW work*/
#define my_type float /**< Input data type*/
#define other_type double/**< The other of double and float*/
#define my_sdf_type SDF_DATATYPE_REAL4 /**< SDF type matching \ref my_type*/
#define MPI_MYTYPE MPI_FLOAT /**< MPI type matching \ref my_type*/
const my_type tiny_my_type=1e-30; /**< Value below which we assume 0*/

#else

#define ADD_FFTW(x) fftw_ ## x /**< Add the correct FFTW library prefix to function/variable names */
#define cplx_type ADD_FFTW(complex)/**< Suitable complex type for FFTW work*/
#define my_type double /**< Input data type*/
#define other_type float/**< The other of double and float*/
#define my_sdf_type SDF_DATATYPE_REAL8 /**< SDF type matching \ref my_type*/
#define MPI_MYTYPE MPI_DOUBLE /**< MPI type matching \ref my_type*/
const my_type tiny_my_type=1e-60; /**< Value below which we assume 0*/

#endif

#define calc_type double /**< Type to do calculations in. This is independent of the input data type \ref my_type */
#define MPI_CALCTYPE MPI_DOUBLE /**< MPI type matching \ref calc_type */

#define tiny_calc_type 1e-12/**< Tiny value for \ref calc_type */

/** @} */
//----------- END TYPE HANDLING----------------------------

/** \defgroup const Constants
*@{ 
\brief Code and physical constants
*
*The subset of these which might change should be logged by each main program on each run using log_code_constants(std::string). They have to be added manually. 
*/
//----------- CONSTANTS ---------------------------------
const size_t MAX_SIZE = 100000;/**< Maximum per-dim array size allowed (per processor if MPI in use) */
const size_t MAX_SIZE_TOT = MAX_SIZE*MAX_SIZE;/**< Maximum overall array size allowed (per processor if MPI in use) */
const int MAX_FILENAME_DIGITS = 15;/**< Maximum number of digits in filename dump number string*/

const int GIT_VERSION_SIZE = 15;/**<Length of git version string*/
const my_type io_verify = 3.0/32.0;/**< An exactly binary representable my_type to verify we're reading what we're writing.*/
const calc_type pi = 3.14159265359;/**< Pi */
const calc_type v0 = 2.997924e8; /**< Speed of light in m/s^2 */
const calc_type me = 9.10938291e-31; /**< Electron mass in kg */
const calc_type mp = me*1836.15267; /**< Proton mass in kg */
const calc_type q0 = 1.602176565e-19; /**< Electron charge in C*/
const calc_type eps0 = 8.85418782e-12; /**< Epsilon_0 permittivity of free space in F/m*/
const calc_type kb = 1.3806488e-23; /**< Boltzman constant J/K*/
const calc_type R_E = 6.371e6; /**<Average Earth radius in m*/

const calc_type GEN_PRECISION = 1e-6;/**< General precision for equality etc*/

const int DEFAULT_N_ANG = 100;/**< Default number of wave normal angles to consider*/
const int ID_SIZE = 10;/**< Length of block ids*/

const int WAVE_WHISTLER = 1; /**< Code to id wave as whistler mode */
const int WAVE_PLASMA = 2; /**< Code to id wave as plasma/Langmuir wave*/
const int WAVE_O = 3; /**< Code to id wave as ordinary EM mode */
const int WAVE_X_UP = 4; /**< Code to id wave as X EM mode, upper branch */
const int WAVE_X_LOW = 5; /**< Code to id wave as X EM mode, lower branch */

const int FUNCTION_NULL = 0; /**< Code to id spectral angular distribution as absent */
const int FUNCTION_DELTA = 1; /**< Code to id spectral angular distribution as delta function (with integral 1) */
const int FUNCTION_GAUSS = 2; /**< Code to id spectral angular distribution as gaussian (with integral 1) */
const int FUNCTION_ISO = 3; /**< Code to id spectral angular distribution as isotropic (with integral 1) */

const std::string OMEGA_CE = "wCe";/**< String specifying  omega_ce in deck.status*/
const std::string OMEGA_PE = "wpe";/**< String specifying  omega_pe in deck.status*/
const std::string DENS_RAT = "dens_rat";/**< String specifying  density ratio in deck.status*/
const std::string DENS_RATH = "dens_rath";/**< String specifying  density ratio in deck.status*/
const std::string DENS = "dens";/**< String specifying  density in deck.status*/
const std::string PPC = "ppc";/**< String specifying ppc in deck.status*/
const std::string VPAR = "vtherm_par";/**<String specifying parallel thermal velocity in deck.status*/
const std::string VPERP = "vtherm_perp";/**<String specifying perpendicular thermal velocity in deck.status*/

const std::string CONSTANTS = " Constant block values after";/**< String denoting start of constant value dump in deck.status*/
const std::string CONSTANTS_END = "Deck state:";/**< String denoting end of constant value dump in deck.status*/

const std::string halp_file = "./files/help/help.txt";/**< Name of command line options help file*/

const std::string LOCAL = "loc";/**< Tag identifying diffusion coefficient processing level: local to space block*/
const std::string BOUNCE_AV = "bav";/**< Tag identifying diffusion coefficient processing level: bounce averaged*/
const std::string GLOBAL = "glb";/**< Tag identifying diffusion coefficient processing level: reduced over all space*/

const char HANDLED_ARG[2] = "*";/**<Dummy string to flag command line arguments as having been handled*/

const my_type V_MIN = 0;/**< Minimum particle velocity for D*/
const my_type V_MAX = 0.99*v0;/**< Maximum particle velocity for D*/
const my_type TAN_MIN = 0.0;/**< Minimum angle (tan theta) for D, spectra etc. Generally should be 0 or -TAN_MAX*/
const my_type TAN_MAX = 4.0;/**< Maximum angle (tan theta) for D, spectra etc*/
const my_type ANG_MIN = 0.0;/**< Minimum angle for spectra etc. Generally should be 0 or -ANG_MAX*/
const my_type ANG_MAX = pi/2;/**< Maximum angle for spectra etc*/
const my_type DEFAULT_SPECTRUM_ANG_STDDEV = 0.2;/**< "Std Dev" of angular distribution (in tan theta) of spectrum*/
const my_type SPECTRUM_THRESHOLD = 1e-3;/**< Fraction of peak power considered to be "significant" spectral power*/

/** @} */
//----------- END CONSTANTS ---------------------------------

/** \defgroup str Data Structures
\brief Data structures for constants etc
*@{ */
//----------- STRUCTURES ---------------------------------

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
  float dens_factor;/**<Ratio of total plasma density to background plasma density*/
};
extern deck_constants my_const;

/** \brief MPI information
*
*Holds info on MPI: processor ranks etc
*/
struct mpi_info_struc{

  int rank;/**< Rank of current processor*/
  int n_procs;/**<Global number of processors*/
};

const struct mpi_info_struc mpi_info_null = {0, -1};/**< Null MPI struct for single threaded jobs, without having to compile the SDF libraries seperately*/

extern mpi_info_struc mpi_info;

/** \brief Reduced refractive index
*
*Contains refractive index mu, the two derivative needed for diffusion coefficient calculation, the phi function of e.g. Albert \cite Albert2005 and error flag
*/
struct mu_dmudom{
  calc_type mu; /**< Refractive index */
  calc_type dmudom; /**< d mu / d omega (wave frequency)*/
  calc_type dmudtheta; /**< d mu / d theta (wave normal angle) */
  calc_type phi;/**< Phi from Albert \cite Albert2005 */
  int err; /**< 0 if mu found successfully, 1 else*/
  calc_type cone_ang;/**<Resonance cone angle*/
};

/** \brief General command line arguments
*
* Processed command line arguments used across several programs
*/
struct setup_args{
  size_t time[3];/**< Start and end dump numbers*/
  bool use_row_time;/**<Whether to use time[2] for sizing*/
  int space[2];/**< Local space block start and end*/
  std::string block;/**< Block ID to use (ex, bz etc)*/
  std::string file_prefix;/**< Prefix part of file names*/
  int n_space;/**< Number of space blocks in global x direction*/
  size_t per_proc;/**< Resulting number of space blocks per proc*/
};

/** \brief Command line arguments for spectra
*
* Processed command line arguments for spectra used across several programs
*/
struct spect_args{
  int fuzz;/**<Fuzz for spectral cutout*/
  int smth;/**<Smoothing width for output spectrum*/
  size_t n_ang;/**<Number of angles for output spectrum*/
  int wave;/**<Wave type ID (see support.h WAVE_* )*/
  int ang;/**< Angular function type (can be FUNCTION_NULL) */
  float ang_sd;/**< Width for angular function (if applicable)*/
  bool mask;/**< Flag to output spectrum extraction mask to file also*/
};

/** \brief D coefficient report
*
*Contains information on D calculation such as resonances used
*/
struct d_report{
  bool error;/**<Whether IO or setup errors occured*/
  size_t n_av;/**< Average n_max used for calc*/
  size_t n_max;/**< Max n_max used in calcs*/
  size_t n_min;/**< Min n_min used in calcs. Not '-' is omitted*/
  bool single_n;/**<Flag showing that a single resonance, n_av, was used*/
};
/** @} */
//-----------END STRUCTURES ---------------------------------
/** @} */

//----------- SPECIFIC HELPER FUNCTIONS -------------
/** \defgroup halp Helper functions
*\brief Groups of helpers of various purpose
*@{*/

/** \defgroup main Main Helper Functions
*@{ 
* \brief General helpers
*
*Contains MPI helpers, argument processing, and some general data handling functions
*/

/********MPI and code helpers ****/
int local_MPI_setup(int argc, char *argv[]);
void safe_exit();
void share_consts();

/********IO helpers ****/
void get_deck_constants(std::string file_prefix);
setup_args process_command_line(int argc, char *argv[]);
spect_args spect_process_command_line(int argc, char *argv[]);
void process_command_line_help_arg(int argc, char *argv[], char help_id);
void print_help(char code=0);
void log_code_constants(std::string file_prefix);
int extract_num_time_part(std::string name);
std::pair<int, int> extract_space_part(std::string name);
std::vector<std::string> process_filelist(int argc, char *argv[]);

class data_array;
/********Data helpers ****/
void divide_domain(std::vector<size_t> dims, size_t space[2], int per_proc, int block_num);
my_type get_ref_Bx(std::string file_prefix, size_t space_in[2], size_t time_0);
data_array get_Bx(std::string file_prefix, size_t space_in[2], size_t time_0);
bool flatten_fortran_slice(my_type * src_ptr, my_type* dest_ptr, size_t n_dims_in, size_t * dims_in, size_t flatten_on_dim,size_t flat_start=0, size_t flat_stop=-1);//I know -1 will overflow, that is what I want

int where(const my_type * ax_ptr, int len, my_type target);

/** Relativistic gamma from velocity
@param v Particle velocity
@return Corresponding relativistic gamma
*/
inline calc_type gamma_rel(calc_type v){
  return 1.0/std::sqrt(1.0 - v*v/v0/v0);
}

/** @} */

//----------- HELPER TYPE FUNCTIONS -------------

/** \defgroup help Global Helper and Maths Functions
*@{ 
*\brief String handling, IO and maths helpers
*/

/********IO helpers ****/
std::string read_wipps_version_string(std::string filename);
bool check_wipps_version(std::string filename);

void my_print(std::string text, int rank, int rank_to_write=0, bool noreturn=false);
void my_print(std::fstream * handle, std::string text, int rank, int rank_to_write=0, bool noreturn=false);
void my_error_print(std::string text, int rank, int rank_to_write=0, bool noreturn=false);
void my_error_print(std::fstream * handle, std::string text, int rank, int rank_to_write=0, bool noreturn=false);

/********String handling helpers ****/
std::string mk_str(int i);/**<Converts int to string*/
std::string mk_str(size_t i); /**< Long int to string*/
std::string mk_str(bool b);/**<Converts bool to string*/
//std::string mk_str(size_t i){ return mk_str((int) i);} /**<Converts size_t to string*/
std::string mk_str(double i, bool noexp=0);/**<Converts double to string*/
std::string mk_str(float i, bool noexp=0);/**<Converts float to string*/
std::string mk_str(long double i, bool noexp=0);/**<Converts long double to string*/
std::string mk_str(char * str);/**<Convert C string to std::string*/
void trim_string(std::string &str, char ch=' '); /**< Trim all leading/trailing ch's from str*/
long checked_strtol(const char * str, bool quiet=false);
float checked_strtof(const char * str, bool quiet=false);
std::string replace_char(std::string str_in, char ch, char repl);/**<Replace all occurences of character ch in string*/
std::string append_into_string(const std::string &in, const std::string &infix);
bool parse_name_val(std::string in, std::string &name, std::string &val);
inline std::string str_to_upper(std::string str){
/** Convert string to upper case */

  for(size_t i=0; i<str.size(); i++){
    if(str[i] >='a' and str[i]<='z') str[i] -=32;
  }
  return str;
}
inline std::string str_to_lower(std::string str){
/** Convert string to lower case */

  for(size_t i=0; i<str.size(); i++){
    if(str[i] >='A' and str[i]<='Z') str[i] +=32;
  }
  return str;
}

int compare_as_version_string(std::string str, std::string vers_str=VERSION, bool minor=false);

/********Maths helpers ****/
template<typename T> T integrator(T * start, int len, T * increment);
template<typename T> void inplace_boxcar_smooth(T * start, int len, int width, bool periodic = 0);
calc_type square_integrator(calc_type * start, int len, calc_type * increment);
std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn);
template<typename T> T interpolate_linear(T axis[2], T vals[2], T target);
template<typename T> T interpolate_nearest(T axis[2], T vals[2], T target);


/********Arithmetic operations ****/
/** \defgroup aux_fns Auxilliary functions
*\brief Named free functions
*
*These can be used with the various array apply functions to do arithmetic on arrays without having to create lambdas
*@{ */

inline my_type subtract(my_type lhs, my_type rhs){return lhs - rhs;}/**< Element-wise subtraction*/
inline my_type add(my_type lhs, my_type rhs){return lhs + rhs;}/**< Element-wise addition*/
inline my_type divide(my_type lhs, my_type rhs){return lhs/rhs;}/**< Element-wise division*/
inline my_type multiply(my_type lhs, my_type rhs){return lhs*rhs;}/**< Element-wise multiplication*/
/** @} */

/** @} */
/** @} */
//----------- END HELPER TYPE FUNCTION DECLARATIONS -----------

#endif
