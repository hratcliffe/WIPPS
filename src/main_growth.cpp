/** \file main_growth.cpp \brief Helper program to get expected and/or actual growth rates of waves
*
* Can take sdf files or derived spectrum files, extract the approximate wave growth rates (peak or integrated) and output these, plus a theoretical rate, OR just the theoretical rate. Parameters are obtained as in Main, so requires a plasma.conf file and deck.status file.
* Depends on the SDF file libraries, the FFTW library. A set of test arguments is supplied. Call using ./growth `<growth_test_pars` to use these.
  \author Heather Ratcliffe \date 11/02/2016.
*/


#include <math.h>
#include <cmath>
#include <fstream>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <iostream>
#include <stdio.h>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>
//FFTW3 Fourier transform libs

#include "support.h"
#include "main_support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "non_thermal.h"


deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/

struct g_args{

  int src;/**<Whether to process no files (analytic only), sdf or spectrum 0, 1,2, respectively*/
  std::string spect_file;/**<File listing spectra if using this option*/
};


calc_type * make_momentum_axis(int n_momenta, calc_type v_max_over_c);

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis)

int main(int argc, char *argv[]){
/**
*In theory, these classes and functions should be named well enough that the function here is largely clear. Remains to be seen, eh?
*
*/

  int err;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  MPI_Barrier(MPI_COMM_WORLD);

  setup_args cmd_line_args = process_command_line(argc, argv);
  g_args extra_args = g_command_line(argc, argv);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/




  if(extra_args.src == 1){
    //Read SDF files and derive...
    my_print("Processing sdf files", mpi_info.rank);
  
  
  
  }else if(extra_args.src == 2){
    //Read spectrum files and derive...
    my_print("Processing spectrum files", mpi_info.rank);
  
  
  
  
  }
  
  //Now we do the analytics.

  my_print("Calculating growth rates", mpi_info.rank);

  plasma * my_plas = new plasma(my_const.omega_ce * me/std::abs(q0), cmd_line_args.file_prefix);

  non_thermal * my_elec = new non_thermal(cmd_line_args.file_prefix);

  const int n_momenta = 1000;
  calc_type * p_axis, *growth_rate;
  
  p_axis = make_momentum_axis(n_momenta, 0.9);
  growth_rate = (calc_type*) malloc(n_momenta, sizeof(calc_type));

  if(!p_axis || ! growth_rate){
    my_print("Cannot allocate memory for arrays", mpi_info.rank);
    
    if(p_axis) free(p_axis);
    if(growth_rate) free(growth_rate);
    //Free whichever was allocated
    
    safe_exit();
  }

  get_growth_rate(my_plas, my_elec, n_momenta, p_axis);


  safe_exit();

}

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis){






}


calc_type * make_momentum_axis(int n_mom, calc_type v_max){

  //Max velocity to consider norm'd to c
  calc_type dp = m0*v_max * v0/ std::sqrt(1.0- v_max*v_max) / (calc_type) n_mom;
  //Momentum step size, including gamma

  calc_type * p_axis;
  p_axis = (calc_type*) malloc(n_mom, sizeof(calc_type));
  if(!p_axis) return nullptr;

  p_axis[0] = 0.0;
  for(int i=1; i< n_momenta; ++i) p_axis[i] = p_axis[i-1] + dp;
  //Set momenta from 0 to v_max. Do this once.

  return p_axis;

}

g_args g_command_line(int argc, char * argv[]){
/** Check whether to process no files (analytic only), sdf or spectrum 0, 1,2, respectively. Looping through again is silly, but we're stealing from main in chunks here... But if we have a -f arg we use sdf files, if a -s we use the spectra listed in that file, if neither, we output analytic only and if both the last one is used*/

  g_args extra_cmd_line;

  extra_cmd_line.src = 0;
  extra_cmd_line.spect_file = "";
  
  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-s")==0 && i < argc-1){
      extra_cmd_line.spect_file = atoi(argv[i+1]);
      extra_cmd_line.src = 2;
    }
    if(strcmp(argv[i], "-f")==0 && i < argc-1) extra_cmd_line.src = 1;
  }
  return extra_cmd_line;


}


