
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

#include "support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "tests.h"

/** \defgroup utils Utility programs
*@{ */

/** \defgroup diff_util Diffusion calculation utility
*@{
*\brief Utility to calculate a particle diffusion coefficient
*
*Calculates a particle diffusion coefficient from given data, in the form of sdf files, ffts or spectrum files. The latter can be created by the generate_ffts and FFT_to_spectrum utils. The resulting particle diffusion coefficient are calculated using Lyons 1974 a, b, Albert 2005 and such. Note that this makes no sense for E fields!
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions. A set of test arguments is supplied. Call using ./calculate_diffusion `<test_pars` to use these. Or try ./calculate_diffusion -h for argument help
\verbinclude help_i.txt
  \author Heather Ratcliffe \date 17/02/2017
*/

const char PER_UTIL_HELP_ID = 'i';

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

  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  setup_args cmd_line_args = process_command_line(argc, argv);

  std::vector<std::string> file_list;
  if(cmd_line_args.is_list){
    file_list = process_filelist(argc, argv);
    cmd_line_args.per_proc = std::ceil( (float) file_list.size() / (float) mpi_info.n_procs);
  }
  
  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /* Get constants from deck and share to other procs*/

  log_code_constants(cmd_line_args.file_prefix);

  my_print("Processing "+mk_str(cmd_line_args.per_proc)+" blocks per core", mpi_info.rank);

  char block_id[ID_SIZE];
  strcpy(block_id, cmd_line_args.block.c_str());

  reader my_reader = reader(cmd_line_args.file_prefix, block_id);
  if(my_reader.current_block_is_accum()) cmd_line_args.use_row_time = true;
  int n_tims;
  if(!cmd_line_args.use_row_time){
    n_tims = std::max((int) (cmd_line_args.time[1]-cmd_line_args.time[0]), 1);
  }else{
    if(cmd_line_args.time[2] == 0){
      my_print("Please specify number of rows to read", mpi_info.rank);
      exit(1);
    }
    n_tims = cmd_line_args.time[2];
  }
  my_reader.update_ref_filenum(cmd_line_args.time[0]);
  size_t my_space[2]={0, 0};

  size_t n_dims;
  std::vector<size_t> dims;
  err = my_reader.read_dims(n_dims, dims);
  if(err) safe_exit();
  int space_dim = dims[0];
  
  controller contr = controller(cmd_line_args.file_prefix);

  //---------------- Now we loop over blocks per proc-------
  for(size_t block_num = 0; block_num<cmd_line_args.per_proc; block_num++){
    
    data_array dat_fft;
    
    if(!cmd_line_args.is_list){
      /* This replaces any -1 in space input with suitable sizes*/
      divide_domain(dims, my_space, cmd_line_args.per_proc, block_num);
      space_dim = my_space[1]-my_space[0];
      
      MPI_Barrier(MPI_COMM_WORLD);
      //--------------THIS will slightly slow down some cores to match the slowest. But it makes output easier. Consider removing if many blocks

      data_array dat = data_array(space_dim, n_tims);

      if(!dat.is_good()){
        my_error_print("Data array allocation failed. Aborting.", mpi_info.rank);
        return 0;
      }

      err = my_reader.read_data(dat, cmd_line_args.time, my_space);
      if(err == 1) safe_exit();

      if(err == 2) n_tims = dat.get_dims(1);
      //Check if we had to truncate data array...

      dat.B_ref = get_ref_Bx(cmd_line_args.file_prefix, my_space, cmd_line_args.time[0] == 0 ? 1: cmd_line_args.time[0]);
      //Get ref B using specfied file but skip 1st ones as they seem broken
      dat_fft = data_array(space_dim, n_tims);
    
      if(!dat_fft.is_good()){
        my_error_print("Data array allocation failed. Aborting.", mpi_info.rank);
        return 0;
      }

      err = fft_array(dat, dat_fft);
      
      if(mpi_info.rank ==0) MPI_Reduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      else MPI_Reduce(&err, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);
    
    }
    else if(!cmd_line_args.is_spect){
      //We just read the specified file
      std::string filename;
      if(cmd_line_args.per_proc*mpi_info.rank + block_num < file_list.size()){
        filename = file_list[cmd_line_args.per_proc*mpi_info.rank + block_num];
        dat_fft = data_array(filename);
 
      }else{
        continue;
        //We've run out of files, must have been non-integer division
      }
    }
    if(!cmd_line_args.is_spect){
      contr.set_plasma_B0(dat_fft.B_ref);
      contr.add_spectrum(space_dim, DEFAULT_N_ANG, true);

      contr.get_current_spectrum()->generate_spectrum(dat_fft);
      //Now we have some test spectral data we can work with...
    }
    else{
      std::string filename;
      if(cmd_line_args.per_proc*mpi_info.rank + block_num < file_list.size()){
        filename = file_list[cmd_line_args.per_proc*mpi_info.rank + block_num];
        contr.add_spectrum(filename);

      }else{
        continue;
        //We've run out of files, must have been non-integer division
      }
    }
    //Now we have some spectral data, either from sdf files, from an FFt dump or from a spectrum dump
    
    contr.add_d(cmd_line_args.d[0], cmd_line_args.d[1]);
    contr.get_current_d()->calculate();


  }
  //-----------------end of per_proc loop---- Now controller holds one spectrum and d per block
  MPI_Barrier(MPI_COMM_WORLD);
  
  bounce_av_data bounce_dat;
  my_space[0] = 0;
  my_space[1] = dims[0];

  //Set up our bounce-averaging params and do the average
  //Accumulated property is for every file in a directory
  bounce_dat.set_Bx_size(my_space[1]);
  bounce_dat.set_Bx(get_Bx(cmd_line_args.file_prefix, my_space, cmd_line_args.time[0], my_reader.current_block_is_accum()));
  bounce_dat.max_latitude = 90.0;
  bounce_dat.L_shell = 4.5;
  contr.bounce_average(bounce_dat);
  
  contr.save_spectra(cmd_line_args.file_prefix);
  contr.save_D(cmd_line_args.file_prefix);

#ifndef NO_FFT
  ADD_FFTW(cleanup());
#endif
  MPI_Finalize();
  //call these last...


  exit(0);
}
/** @} */
/** @} */


