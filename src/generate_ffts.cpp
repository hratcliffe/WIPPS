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
#include "data_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "main_support.h"

deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/

struct gen_cmd_line{
  int flat_dim;
  bool flat_fft;
  my_type flat_fft_min;
  my_type flat_fft_max;
  std::vector<my_type> limits;
};

gen_cmd_line special_command_line(int argc, char *argv[]);

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

  gen_cmd_line extra_args = special_command_line(argc, argv);
  setup_args cmd_line_args = process_command_line(argc, argv);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/

  my_print("Processing "+mk_str(cmd_line_args.per_proc)+" blocks per core", mpi_info.rank);
  my_print("Input "+cmd_line_args.file_prefix, mpi_info.rank);

  char block_id[ID_SIZE];
  strcpy(block_id, cmd_line_args.block.c_str());

  reader * my_reader = new reader(cmd_line_args.file_prefix, block_id);
  if(my_reader->current_block_is_accum()) cmd_line_args.use_row_time = true;

  int n_tims;
  if(!cmd_line_args.use_row_time){
    n_tims = std::max(cmd_line_args.time[1]-cmd_line_args.time[0], 1);
  }else{
    n_tims = cmd_line_args.time[2];
  }
  my_reader->set_ref_filenum(cmd_line_args.time[0]);

  int my_space[2];
  my_space[0] = cmd_line_args.space[0];
  my_space[1] = cmd_line_args.space[1];

  size_t n_dims;
  std::vector<size_t> dims;
  err = my_reader->read_dims(n_dims, dims);
  if(err){
    my_print("File read error, aborting", mpi_info.rank);
    safe_exit();
  }
  int space_dim = dims[0];
  
  controller * contr;
  contr = new controller(cmd_line_args.file_prefix);

  std::fstream logfile;
  std::string logfilename;
  logfilename = cmd_line_args.file_prefix + "generate_ffts.log";
  logfile.open(logfilename.c_str(),std::ios::app);

  //---------------- Now we loop over blocks per proc-------
  for(int block_num = 0; block_num<cmd_line_args.per_proc; block_num++){

    divide_domain(dims, my_space, cmd_line_args.per_proc, block_num);
    space_dim = my_space[1]-my_space[0];
    
    MPI_Barrier(MPI_COMM_WORLD);
    //--------------THIS will slightly slow down some cores to match the slowest. But it makes output easier. Consider removing if many blocks
    data_array dat;
    if(n_dims ==1 || (n_dims ==2 && my_reader->current_block_is_accum())){
      dat = data_array(space_dim, n_tims);
    }else if(n_dims == 2){
      dat = data_array(space_dim, dims[1], n_tims);
    
    }else if(n_dims == 3){
      dat = data_array(space_dim, dims[1], dims[2], n_tims);
    
    }else{
      my_print("More than 3-D input data not supported by this utility", mpi_info.rank);
    }//Here I have no more than 3-D in space
  
    if(!dat.is_good()){
      my_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = my_reader->read_data(dat, cmd_line_args.time, my_space);
    /** \todo Have this attempt the next file before failing*/
    //IMPORTANT: we use the same array for each space block. So if we had insufficient files etc this will truncate the array ONLY ONCE
    
    if(err == 1){
      my_print("Data read failed. Aborting", mpi_info.rank);
      safe_exit();
    }
    if(err == 2) n_tims = dat.get_dims(1);
    //Check if we had to truncate data array...

    dat.B_ref = get_ref_Bx(cmd_line_args.file_prefix, my_space, cmd_line_args.time[0] == 0 ? 1:cmd_line_args.time[0], my_reader->current_block_is_accum());
    //Get ref B using specfied file but skip 0th ones as they seem broken

    //Checks limits and denormalise freq. if applicable
    std::vector<my_type> lims;
    if(extra_args.limits.size() != 0){
      if(extra_args.limits.size() != 2*dat.get_dims()) my_print("Please supply 2 limits per dimension. Output will be untrimmed");
      else{
        lims = extra_args.limits;
        if(extra_args.flat_dim != dat.get_dims()-1){
          lims[lims.size()-2] *= (dat.B_ref*q0/me);
          lims[lims.size()-1] *= (dat.B_ref*q0/me);
        }
      }
    }
    if(extra_args.flat_dim>=0 && dat.get_dims()>1 && !extra_args.flat_fft){
      dat = dat.total(extra_args.flat_dim);
    }
    
    data_array dat_fft;
    dat_fft.clone_empty(dat);

    if(!dat_fft.is_good()){
      my_print("Data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = dat.fft_me(dat_fft);

    if(mpi_info.rank ==0) MPI_Reduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(&err, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);
    
    //Construct filename. Since the MPI is using block-wise domain decomposition, different processors can't overlap on blocks
    std::string filename, time_str;
    time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(n_tims);
    std::string block = block_id;
    filename = cmd_line_args.file_prefix+"FFT_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
    std::fstream file;
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()){
      if(lims.size() == 2*dat.get_dims()){
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
  
  //contr->save_spectra(cmd_line_args.file_prefix);

  logfile.close();
  //Cleanup objects etc
  delete my_reader;
  delete contr;

  std::cout<<"Grep for FAKENUMBERS !!!!"<<std::endl;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...


  exit(0);
}

gen_cmd_line special_command_line(int argc, char *argv[]){
//Do special command line processing here

  gen_cmd_line values;
  values.flat_dim = -1;
  values.flat_fft = false;
  values.flat_fft_min = 0.0;
  values.flat_fft_max = 0.0;
  
  for(int i=0; i< argc; i++){
    if(strcmp(argv[i], "-h")==0){
      print_help();
      print_help('g');
      exit(0);
    }
    else if(strcmp(argv[i], "-flat_dat")==0 && i < argc-1) values.flat_dim = atoi(argv[i+1]);
    else if(strcmp(argv[i], "-flat_fft")==0  && i < argc-3){
      values.flat_dim = atoi(argv[i+1]);
      values.flat_fft = true;
      values.flat_fft_min = atof(argv[i+2]);
      values.flat_fft_max = atof(argv[i+3]);
    }
    else if(strcmp(argv[i], "-lims")==0 && i < argc-1){
      while(i<argc-1 && (argv[i+1][0]!= '-'  || ((argv[i+1][1] >='0' && argv[i+1][1] <='9') || argv[i+1][1] =='.'))){
        //Checks if next argument is a new flag, but allows negative numbers
        values.limits.push_back(atof(argv[i+1]));
        i++;
      }
    }
  }
  
  return values;
}
