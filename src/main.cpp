/** \file main.cpp \brief Main program
*
* Calculates a diffusion coefficient from data. Data can be either a range of SDF files or a list of FFT or spectrum files, from e.g. supplied generate_ffts utility. The resulting particle diffusion coefficient are calculated using Lyons 1974 a, b, Albert 2005 and such. Note that this makes no sense for E fields!
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions. A set of test arguments is supplied. Call using ./main `<test_pars` to use these. Or try ./main -h for argument help
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
#include "tests.h"

#ifdef RUN_TESTS_AND_EXIT
tests* test_bed;/**<Test bed for testing */
#endif
//We wrap in ifdef for nice Doxygen docs

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

#ifdef RUN_TESTS_AND_EXIT
  my_print("Running basic tests", mpi_info.rank);

  if(mpi_info.rank == 0) get_deck_constants("./files/test");
  share_consts();

  test_bed = new tests();
  test_bed->set_verbosity(2);
  err = test_bed->run_tests();
  delete test_bed;
  exit(err);

#else
  //Explain the code....
  
  my_print("Welcome to ...", mpi_info.rank);
  print_help();
  
#endif

  MPI_Finalize();
  
  exit(0);
}

