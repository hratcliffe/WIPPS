
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

#ifndef NO_FFT
#include <fftw3.h>
//FFTW3 Fourier transform libs
#endif

#include "main.h"
#include "support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"


/** \defgroup utils Available programs
*@{ 
*\brief Top-level programs available
*
*Contains programs to perform FFTs, create spectra, calculate growth etc. Build with `make utils` List with `make list_utils` Command line argument help is available using ./{util_name} -h
*/

/** \defgroup main_prog Info program
*@{
\brief Info program
*
* This prints information about the code or compiled in test mode it runs the tests and prints results.
*/

#ifndef RUN_TESTS_AND_EXIT
const char PER_UTIL_HELP_ID = ' ';/**<ID to identify help file for this utility*/
#endif

/** \brief Main program
*
* Compiled in test MODE this runs the testing code. Otherwise it prints info.
  \author Heather Ratcliffe \date 18/09/2015.
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char *argv[]){
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);
    
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef RUN_TESTS_AND_EXIT
  my_print("Running tests", mpi_info.rank);

  if(mpi_info.rank == 0) get_deck_constants("./files/test");
  share_consts();

  test_bed = new tests();
  test_bed->set_verbosity(2);
  test_bed->set_runtime_flags(argc, argv);
  int err = test_bed->run_tests();
  delete test_bed;
  exit(err);

#else
  //Print general code help
  print_help(PER_UTIL_HELP_ID);
  
#endif

  MPI_Finalize();
  
  exit(0);
}
/**@}*/
/**@}*/
