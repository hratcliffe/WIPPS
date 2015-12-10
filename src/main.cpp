/** \file main.cpp \brief Main program
*
* This should: open a sequence of SDF files using the SDF library and read E and B field data. Fourier transform it. Extract frequency/wavenumber and angular spectra (if possible). Calculate the resulting particle diffusion coefficients using Lyons 1974 a, b, Albert 2005 and such.
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions.
  \author Heather Ratcliffe \date 18/09/2015.
*/


#include <math.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>
//FFTW3 Fourier transform libs

#include "main.h"
#include "support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "tests.h"

using namespace std;

deck_constants my_const;/**< Physical constants*/
mpi_info_struc mpi_info;/**< MPI data */

tests* test_bed;/**Test bed for testing */

void get_deck_constants();

int main(int argc, char *argv[]){
/**
*In theory, these classes and functions should be named well enough that the function here is largely clear. Remains to be seen, eh?
*
*/

  int ierr,err;
  
{
  int rank, n_procs;

  ierr = MPI_Init(&argc, &argv);
  //Note any other command line arg processing should account for this...

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

  mpi_info.rank = rank;
  mpi_info.n_procs = n_procs;
}




  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  get_deck_constants();
  /** Get constants from deck*/

#ifdef RUN_TESTS_AND_EXIT
  my_print("Running basic tests", mpi_info.rank);
  test_bed = new tests();
  test_bed->set_verbosity(2);
  test_bed->run_tests();
  delete test_bed;

  return 0;
#else

  //Actually do the code...

  char block_id[10]= "ex";
  reader * my_reader = new reader("", block_id);

  int tim_in[2], space_in[2];
  tim_in[0]=0;
  tim_in[1]=10;
  space_in[0]=0;
  space_in[1]=-1;
  /** \todo These should be specifiable via command line args. Or an input file. Or in parallle we want to subdivide into them
*/

  int n_tims = max(tim_in[1]-tim_in[0], 1);

  int n_dims;
  std::vector<int> dims;
  my_reader->read_dims(n_dims, dims);

  if(n_dims !=1) return 1;
  /**for now abort if data file wrong size... \todo FIX*/

  data_array  * dat = new data_array(dims[0], n_tims);
  data_array * dat_fft = new data_array(dims[0], n_tims);

  if(!dat->is_good() or !dat_fft->is_good()){
    my_print("Bugger, data array allocation failed. Aborting.", mpi_info.rank);
    return 0;
  }

  my_reader->read_data(dat, tim_in, space_in);

  err = dat->fft_me(dat_fft);
  my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);

  fstream file;
  file.open("Tmp.txt", ios::out|ios::binary);
  if(file.is_open()) dat->write_to_file(file);
  file.close();

  controller * contr;
  contr = new controller();

  int row_lengths[2];
  row_lengths[0] = 4096;
  row_lengths[1] = DEFAULT_N_ANG;
  
  contr->add_spectrum(row_lengths, 2);
  contr->get_current_spectrum()->make_test_spectrum();

  file.open("Tmp_spectrum.txt", ios::out|ios::binary);

  if(file.is_open() && contr->get_current_spectrum()) contr->get_current_spectrum()->write_to_file(file);
  file.close();

  //Now we have some test spectral data we can work with...

  contr->add_d(10, 10);
  contr->get_current_d()->calculate();

  //Cleanup objects etc
  delete my_reader;
  delete dat;
  delete dat_fft;
  delete contr;

  cout<<"Grep for FAKENUMBERS !!!!"<<endl;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...
#endif

  exit(0);
}


int where(my_type * ax_ptr, int len, my_type target){
/** \todo handle sign etc in here...
*/
  int whereb(my_type * ax_ptr, int len, my_type target, int &cut,int sign); //Recursive function to do the finding

  int sign = 1.0;
  int cut = 0;
  if(ax_ptr[0] >=target || ax_ptr[len-1] < target) return -1;
  return whereb(ax_ptr, len, target, cut, sign);

}

int whereb(my_type * ax_ptr, int len, my_type target,int &cut, int sign){
/**\brief Recursive binary bisection find. 
*
*First index where ax_ptr exceeds target. Note scoped only to where()
*/
  
  if(len==1){
     cut++;
     return cut;
  }
  else if(sign*ax_ptr[len/2] >= sign*target){
    whereb(ax_ptr, len/2, target, cut, sign);

    return cut;
  }
  else if(sign*ax_ptr[len/2] < sign*target){
    cut+= len/2;
    whereb(ax_ptr+len/2, len-len/2, target, cut, sign);
    return cut;
  }

  return -1;
  //Shouldn't ever reach this case, but.
}



void get_deck_constants(){
/** \brief Setup run specific constants
*
*This will read deck.status and parse values for user defined constants etc. It will rely on using the specific deck, because it has to look for names. Any changes to deck may need updating here. \todo Write it. \todo Document H assumption...
*/

my_const.omega_ce = 17588.200878;
my_const.omega_ci = my_const.omega_ce * me/mp;
//assumes same charge magnitude, I.e. H plasma

my_const.omega_pe = 5.0*35176.401757;

}

void my_print(std::string text, int rank, int rank_to_write){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also.
*/
  if(rank == rank_to_write){
  
    std::cout<< text<<std::endl;
  }

}
void my_print(fstream * handle, std::string text, int rank, int rank_to_write){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also.
*/
  if(rank == rank_to_write && handle!=nullptr){
    *handle<<text<<std::endl;
  }else if(rank == rank_to_write){
    std::cout<<text<<std::endl;

  }

}


std::string mk_str(int i){

  char buffer[25];
  std::sprintf(buffer, "%i", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(double i){

  char buffer[25];
  std::snprintf(buffer, 25, "%e", i);
  std::string ret = buffer;
  return ret;
  
}
std::string mk_str(float i){

  char buffer[25];
  std::snprintf(buffer, 25, "%e", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(bool b){

  if(b) return "1";
  else return "0";

}

std::string mk_str(long double i){return mk_str((double) i);};

template<typename T> T integrator(T * start, int len, T * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  T value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template float integrator<float>(float *, int, float *);
template double integrator<double>(double *, int, double *);
//We need both float and double versions

calc_type square_integrator(calc_type * start, int len, calc_type * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  calc_type value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i]*start[i] + start[i+1]*start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template<typename T> void inplace_boxcar_smooth(T * start, int len, int width, bool periodic){
/** \brief Boxcar smoothing of specified width
*
*Smooths the array given by start and len using specified width. If periodic is set the ends wrap around. Otherwise they one-side
*/

  if(width > len) my_print("Really? How am I meant to smooth that?", mpi_info.rank);
  calc_type result = 0.0;
  int edge = width/2;
  for(int i=0; i<width; ++i) result += start[i];
  for(int i=edge+1; i<len-edge; ++i){
    result -= start[i-edge-1];
    result += start[i+edge];
    //remove behind and add in front. Faster for width>2
    start[i] = result / (calc_type) width;
  }

  //Handle ends
  if(periodic){
    result = 0.0;
    for(int i=0; i<width-1; ++i) result += start[i];
    //first width-1
    int wrap = len -1;
    //and one wrapped back
    result +=start[wrap];
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - i - 1];
      start[edge-i] = result / (calc_type) width;
    }
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - edge - i - 1];
      start[len-i-1] = result / (calc_type) width;
    }


  }else{
    //we just truncate at the bottom
    
  
  
  }




}
template void inplace_boxcar_smooth(calc_type *, int, int, bool);

std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn){
/** \brief Finds roots of cubic x^3 + an x^2 + bn x + cn = 0
*
* Uses Num. Rec. equations, which are optimised for precision. Note that if x >>1 precision errors may result. Returns real solutions only
*/

  calc_type Q, R, bigA, bigB, Q3, R2, bigTheta;
  std::vector<calc_type> ret_vec;

  Q = (std::pow(an, 2) - 3.0 * bn)/9.0;
  R = (2.0* std::pow(an, 3) - 9.0 * an *bn + 27.0*cn)/54.0;
  
  R2 = std::pow(R, 2);
  Q3 = std::pow(Q, 3);
  
  if( R2 < Q3){
    
    bigTheta = std::acos(R/sqrt(Q3));
    calc_type minus2sqrtQ = -2.0*std::sqrt(Q);
    
    ret_vec.push_back(minus2sqrtQ*std::cos(bigTheta/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta + 2.0*pi)/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta - 2.0*pi)/3.0) - an/3.0);

  }else{
    calc_type ret_root;
    bigA = - boost::math::sign(R)*std::pow((std::abs(R) + std::sqrt(R2 - Q3)), 1.0/3.0 );

    (bigA != 0.0) ? (bigB = Q / bigA) : (bigB = 0.0);
    ret_root = (bigA + bigB) - an/3.0;

    ret_vec.push_back(ret_root);
  }

/** Used to test when writing
  calc_type tmp;
  for(int i=0; i<ret_vec.size(); ++i){
    
    tmp = std::pow(ret_vec[i], 3) + an*std::pow(ret_vec[i], 2) + bn*ret_vec[i] + cn;
    std::cout<<"solution gives "<<tmp<<std::endl;
  
  }
*/
  return ret_vec;

}

template<typename T> T interpolate(T* axis, T* vals, T target, int pts){
/** Interpolate vals on axis to target value
*
*For pts=1 uses closest value, pts=2 uses 2 pt linear, \todo add more pts options
*/

  T ret = 0.0;
  if(pts ==1){
    //select closer value
    if(std::abs(target - axis[0]) <= std::abs(target - axis[1])) ret = vals[0];
    else ret = vals[1];
  
  }else if(pts ==2){
    ret = (std::abs(target - axis[1]) * vals[0] + std::abs(target - axis[0]) * vals[1])/(std::abs(axis[1] - axis[0]));
  
  }

  return ret;
}
template float interpolate(float*, float*, float, int);
template double interpolate(double*, double*, double, int);

