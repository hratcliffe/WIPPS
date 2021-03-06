
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
#include "data_array.h"
#include "d_coeff.h"
#include "spectrum.h"

/** \ingroup utils */

/** \defgroup fft_util FFT generator utility
*@{ 
*\brief Utility to read files and perform Fourier transforms
*
*Opens files, extracts specified fields, does FFT, trims to specified boundaries and writes to file. Can optionally flatten the raw data before Ft-ing or the FT-d data before output. This routine can use multiple cores to process seperate spatial blocks.
\verbinclude help_g.txt
  \author Heather Ratcliffe \date 04/07/2016

*/

const char PER_UTIL_HELP_ID = 'g';/**<ID to identify help file for this utility*/

struct gen_cmd_line{
  int flat_dim;/**<Flattening dimension number*/
  bool flat_fft;/**<Flatten after FFT*/
  my_type flat_fft_min;/**<Lower band limit for FFT flattening*/
  my_type flat_fft_max;/**<Upper band limit for FFT flattening*/
  std::vector<my_type> limits;/**< Limits to trim output FFT to*/
};/**<\brief Additional command line arguments for FFT generation utility*/

gen_cmd_line special_command_line(int argc, char *argv[]);

/** \brief Main program
*
* Generate an FFT
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */
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

  //Get the special args and then the rest
  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  gen_cmd_line extra_args = special_command_line(argc, argv);
  setup_args cmd_line_args = process_command_line(argc, argv);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /* Get constants from deck and share to other procs*/

  log_code_constants(cmd_line_args.file_prefix);
  /* Do this so we have the spectrum wavetype and angle codes*/

  my_print("Processing "+mk_str(cmd_line_args.per_proc)+" blocks per core", mpi_info.rank);
  my_print("Input "+cmd_line_args.file_prefix, mpi_info.rank);

  reader my_reader = reader(cmd_line_args.file_prefix, cmd_line_args.block,cmd_line_args.time[0]);
  if(my_reader.current_block_is_accum()) cmd_line_args.use_row_time = true;

  //Get all dimensions
  int n_tims;
  if(!cmd_line_args.use_row_time){
    n_tims = std::max((int)(cmd_line_args.time[1]-cmd_line_args.time[0]), 1);
  }else{
    n_tims = cmd_line_args.time[2];
  }
  size_t my_space[2];
  my_space[0] = cmd_line_args.space[0];
  my_space[1] = cmd_line_args.space[1];

  size_t n_dims;
  std::vector<size_t> dims;
  err = my_reader.read_dims(n_dims, dims);
  if(err){
    my_error_print("File read error, aborting", mpi_info.rank);
    safe_exit();
  }
  int space_dim = dims[0];
  
  controller * contr;
  contr = new controller(cmd_line_args.file_prefix);

  //Open output logfile
  std::fstream logfile;
  std::string logfilename;
  logfilename = cmd_line_args.file_prefix + "generate_ffts.log";
  logfile.open(logfilename.c_str(),std::ios::app);

  //---------------- Now we loop over blocks per proc-------
  for(size_t block_num = 0; block_num<cmd_line_args.per_proc; block_num++){

    divide_domain(dims, my_space, cmd_line_args.per_proc, block_num);
    space_dim = my_space[1]-my_space[0];
    
    MPI_Barrier(MPI_COMM_WORLD);
    //--------------THIS will slightly slow down some cores to match the slowest. But it makes output easier. Consider removing if many blocks
    data_array dat;
    if(n_dims ==1 || (n_dims ==2 && my_reader.current_block_is_accum())){
      dat = data_array(space_dim, n_tims);
    }else if(n_dims == 2){
      dat = data_array(space_dim, dims[1], n_tims);
    
    }else if(n_dims == 3){
      dat = data_array(space_dim, dims[1], dims[2], n_tims);
    
    }else{
      my_error_print("More than 3-D input data not supported by this utility", mpi_info.rank);
    }//Here I have no more than 3-D in space
  
    if(!dat.is_good()){
      my_error_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = my_reader.read_data(dat, cmd_line_args.time, my_space);
    //IMPORTANT: we use the same array for each space block. So if we had insufficient files etc this will truncate the array ONLY ONCE
    
    if(err == 1){
      my_error_print("Data read failed. Aborting", mpi_info.rank);
      safe_exit();
    }
    if(err == 2) n_tims = dat.get_dims(n_dims);
    //Check if we had to truncate data array...

    dat.B_ref = get_ref_Bx(cmd_line_args.file_prefix, my_space, cmd_line_args.time[0] == 0 ? 1:cmd_line_args.time[0]);
    //Get ref B using specified file but skip 0th dumps as they seem broken

    //Checks limits and denormalise freq. if applicable
    std::vector<my_type> lims;
    if(extra_args.limits.size() != 0){
      int do_flat = (extra_args.flat_dim>=0);
      if(extra_args.limits.size() != 2*(dat.get_dims() - do_flat)) my_print("Please supply 2 limits per dimension. Output will be untrimmed", mpi_info.rank);
      else{
        lims = extra_args.limits;
        //Check for -1 to avoid overflow in second comparison
        if(extra_args.flat_dim == -1 || (extra_args.flat_dim != -1 && (size_t) extra_args.flat_dim != dat.get_dims()-1)){
          lims[lims.size()-2] *= (dat.B_ref*q0/me);
          lims[lims.size()-1] *= (dat.B_ref*q0/me);
        }
      }
    }
    //Flatten data pre-fft if requested
    if(extra_args.flat_dim >= 0 && dat.get_dims() > 1 && !extra_args.flat_fft){
      dat = dat.average(extra_args.flat_dim);
    }
    //Setup the fft output array
    data_array dat_fft;
    dat_fft.clone_empty(dat);
    if(!dat_fft.is_good()){
      my_error_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = fft_array(dat, dat_fft);

    if(mpi_info.rank ==0) MPI_Reduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(&err, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);
    
    
    if(extra_args.flat_fft && extra_args.flat_dim >=0){
      //Flatten fft between limits if requested
      dat_fft = dat_fft.total(extra_args.flat_dim, extra_args.flat_fft_min, extra_args.flat_fft_max);
    
    }
    
    //Construct filename. Since the MPI is using block-wise domain decomposition, different processors can't overlap on blocks
    std::string filename, time_str, flat_tag = "";
    if(extra_args.flat_dim >=0 && !extra_args.flat_fft) flat_tag ="_f"+ mk_str(extra_args.flat_dim);
    else if(extra_args.flat_dim >=0) flat_tag ="_fFFT"+ mk_str(extra_args.flat_dim);
    
    time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(n_tims);
    filename = cmd_line_args.file_prefix+"FFT_"+ cmd_line_args.block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + flat_tag+".dat";
    std::fstream file;
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()){
      if(lims.size() == 2*dat_fft.get_dims()){
        dat_fft.write_section_to_file(file, lims);
      }else{
        dat_fft.write_to_file(file);
      }
    }
    file.close();
    my_print( "FFT section output in "+filename, mpi_info.rank);
    if(logfile) my_print(&logfile, "FFT section output in "+filename, mpi_info.rank);

  }
  //-----------------end of per_proc loop---- Now controller holds one spectrum and d per block
  MPI_Barrier(MPI_COMM_WORLD);
  
  logfile.close();
  //Cleanup objects etc
  delete contr;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...

  exit(0);
}

gen_cmd_line special_command_line(int argc, char *argv[]){
/** \brief Process special arguments to generate_ffts
*
*Part of argument handling is shared with calculate_diffusion so we handle only the extras here and must pass the rest on. So we nullify those we handle here to enable warning for unknown arguments by setting the first character of them all to HANDLED_ARG
*/

  gen_cmd_line values;
  //Default values
  values.flat_dim = -1;
  values.flat_fft = false;
  values.flat_fft_min = 0.0;
  values.flat_fft_max = 0.0;
  
  int i = 1;
  while(i < argc){
    if(strcmp(argv[i], "-flat_dat")==0 && i < argc-1){
      values.flat_dim = checked_strtol(argv[i+1]);
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      i++;
    }
    else if(strcmp(argv[i], "-flat_fft")==0  && i < argc-3){
      values.flat_dim = checked_strtol(argv[i+1]);
      values.flat_fft = true;
      values.flat_fft_min = checked_strtof(argv[i+2]);
      values.flat_fft_max = checked_strtof(argv[i+3]);
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      strcpy(argv[i+2], HANDLED_ARG);
      strcpy(argv[i+3], HANDLED_ARG);
      
      i+=3;
    }
    else if(strcmp(argv[i], "-lims")==0 && i < argc-1){
      while(i<argc-1 && (argv[i+1][0]!= '-'  || ((argv[i+1][1] >='0' && argv[i+1][1] <='9') || argv[i+1][1] =='.'))){
        //Checks if next argument is a new flag, but allows negative numbers
        values.limits.push_back(checked_strtof(argv[i+1]));
        strcpy(argv[i], HANDLED_ARG);
        i++;
      }
      strcpy(argv[i], HANDLED_ARG);

    }
    i++;
  }
  
  return values;
}
/** @} */
/** @} */
