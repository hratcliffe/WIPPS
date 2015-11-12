/** \file main.cpp \brief Main program
*
* This should: open a sequence of SDF files using the SDF library and read E and B field data. Fourier transform it. Extract frequency/wavenumber and angular spectra (if possible). Calculate the resulting particle diffusion coefficients using Lyons 1974 a, b, Albert 2005 and such.
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions.
  \author Heather Ratcliffe \date 18/09/2015.
*/


#include <math.h>
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
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "tests.h"

using namespace std;

deck_constants my_const;/*< Physical constants*/
mpi_info_struc mpi_info;
void get_deck_constants();

void test_bes();
tests* test_bed;


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

#ifdef RUN_TESTS_AND_EXIT
  cout<<"Running basic tests"<<endl;
  test_bed = new tests();
  test_bed->run_tests();
  delete test_bed;

  return 0;
#else

  //Actually do the code...
  get_deck_constants();

  char block_id[10]= "ex";
  reader * my_reader = new reader("", block_id);

  int tim_in[2], space_in[2];
  tim_in[0]=0;
  tim_in[1]=10;
  space_in[0]=0;
  space_in[1]=-1;
  /** \todo These should be specifiable via command line args. Or in parallle we want to subdivide into them
*/

  int n_tims = max(tim_in[1]-tim_in[0], 1);

  int n_dims;
  std::vector<int> dims;
  my_reader->read_dims(n_dims, dims);

  if(n_dims !=1) return 1;
  //for now abort if data file wrong size...

  data_array  * dat = new data_array(dims[0], n_tims);
  data_array * dat_fft = new data_array(dims[0], n_tims);

  if(!dat->data or !dat_fft->data){
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


  /*  file.open("Tmp_dat.txt", ios::in|ios::binary);
    dat.read_from_file(file);
    file.close();
  */

  int row_lengths[2];
  row_lengths[0] = 4096;
  row_lengths[1] = DEFAULT_N_ANG;
  
  spectrum * spect = new spectrum(row_lengths, 2);
  spect->make_test_spectrum();

  file.open("Tmp_spectrum.txt", ios::out|ios::binary);
  if(file.is_open()) spect->write_to_file(file);
  file.close();

  //Now we have some test spectral data we can work with...

// if(sizeof(input_type) != sizeof(calc_type)). Hmm, now we need the data in a higher precision format...
  
  diffusion_coeff * D;
  D = new diffusion_coeff(100, 100);
  plasma * plasm;
  plasm = new plasma();

  D->calculate(spect, plasm);

  //Cleanup objects etc
  delete spect;
  delete my_reader;
  delete dat;
  delete dat_fft;
  delete plasm;
  delete D;

  cout<<"Grep for FAKENUMBERS !!!!"<<endl;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...
#endif
}


void test_bes(){
//test code using bessel funtion. Output values should be zeros.
//TODO expnd these things so we can test what we use before trying to proceed....


//cyl_bessel_j(v, x) = Jv(x)
//cyl_neumann(v, x) = Yv(x) = Nv(x)
cout<<"Bessel test: "<<endl;
double bess, arg;
int index;

index = 0;
arg = 2.40482555769577;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index = 1;
arg=7.01558666981561;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index=5;
arg=12.3386041974669;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index =0;
arg =1.0;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess-0.7651976865579665514497<<endl;

}

void get_deck_constants(){
/** \brief Setup run specific constants
*
*This will read deck.status and parse values for user defined constants etc. It will rely on using the specific deck, because it has to look for names. Any changes to deck may need updating here. \todo Write it. \todo Document H assumption...
*/

my_const.omega_ce = 17588.200878;
my_const.omega_ci = my_const.omega_ce * me/mp;
//assumes same charge magnitude, I.e. H plasma

my_const.omega_pe = 35176.401757;

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
  if(rank == rank_to_write){
    *handle<<text<<std::endl;
  }

}


std::string mk_str(int i){

  char buffer[25];
  std::sprintf(buffer, "%i", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(bool b){

  if(b) return "1";
  else return "0";

}

calc_type integrator(calc_type * start, int len, calc_type * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  calc_type value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }

 return value;

}
