/** \file generate_ffts.cpp \brief Helper program to extract FFTS
*
*
* Opens files, extracts specified fields, does FFTS, trims to specified boundaries and writes to file
  \author Heather Ratcliffe \date 04/07/2016.
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
#include "main_support.h"

deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/


int main(int argc, char *argv[]){

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
  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/

  my_print("Processing "+mk_str(cmd_line_args.per_proc)+" blocks per core", mpi_info.rank);

  char block_id[ID_SIZE];
  strcpy(block_id, cmd_line_args.block.c_str());

  reader * my_reader = new reader(cmd_line_args.file_prefix, block_id);

  int n_tims = std::max(cmd_line_args.time[1]-cmd_line_args.time[0], 1);

  int my_space[2];
  my_space[0] = cmd_line_args.space[0];
  my_space[1] = cmd_line_args.space[1];

  int n_dims;
  std::vector<int> dims;
  err = my_reader->read_dims(n_dims, dims);
  if(err) safe_exit();
  int space_dim = dims[0];
  
  if(n_dims !=1) return 1;
  /**for now abort if data file wrong size... \todo FIX*/

  controller * contr;
  contr = new controller(cmd_line_args.file_prefix);


  //---------------- Now we loop over blocks per proc-------
  for(int block_num = 0; block_num<cmd_line_args.per_proc; block_num++){

    divide_domain(dims, my_space, cmd_line_args.per_proc, block_num);
    space_dim = my_space[1]-my_space[0];
    
    MPI_Barrier(MPI_COMM_WORLD);
    //--------------THIS will slightly slow down some cores to match the slowest. But it makes output easier. Consider removing if many blocks

    data_array  * dat = new data_array(space_dim, n_tims);

    if(!dat->is_good()){
      my_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = my_reader->read_data(dat, cmd_line_args.time, my_space);
    if(err == 1) safe_exit();

    if(err == 2) n_tims = dat->get_dims(1);
    //Check if we had to truncate data array...
    data_array * dat_fft = new data_array(space_dim, n_tims);

    if(!dat_fft->is_good()){
      my_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = dat->fft_me(dat_fft);
    
    if(mpi_info.rank ==0) MPI_Reduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(&err, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);

    int row_lengths[2];
    row_lengths[0] = space_dim;
    row_lengths[1] = DEFAULT_N_ANG;
    
    contr->add_spectrum(row_lengths, 2);
    contr->get_current_spectrum()->make_test_spectrum(cmd_line_args.time, my_space);
    

    std::vector<my_type> lims;
    if(dat.n_dims >=3){
      lims.push_back(-0.002);
      lims.push_back(0.002);
    }
    if(dat.n_dims >=2){
      lims.push_back(-0.002);
      lims.push_back(0.002);
      lims.push_back(-2*my_const.omega_ce);
      lims.push_back(2*my_const.omega_ce);
    
    }
  //Set cutout limits on FFT
    
    std::string block = block_id;
    std::string filename = cmd_line_args.file_prefix+"FFT_"+block +"_"+mk_str(cmd_line_args.time[0])+"_"+mk_str(cmd_line_args.time[1])+"_"+mk_str(cmd_line_args.space[0])+"_"+mk_str(cmd_line_args.space[1]) + ".dat";
    std:::fstream file;
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()) dat_fft->write_section_to_file(file, lims);
    file.close();


    delete dat;

  }
  //-----------------end of per_proc loop---- Now controller holds one spectrum and d per block
  MPI_Barrier(MPI_COMM_WORLD);
  
  contr->save_spectra(cmd_line_args.file_prefix);

  //Cleanup objects etc
  delete my_reader;
  delete contr;

  std::cout<<"Grep for FAKENUMBERS !!!!"<<std::endl;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...
#endif

  exit(0);
}
