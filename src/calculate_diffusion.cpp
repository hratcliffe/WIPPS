
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

#include "support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "tests.h"

/** \ingroup utils */
/** \defgroup diff_util Diffusion calculation utility
*@{
*\brief Utility to calculate a particle diffusion coefficient
*
*Calculates a particle diffusion coefficient from given data, in the form of sdf files, ffts or spectrum files. The latter can be created by the generate_ffts and FFT_to_spectrum utils. Note that cross-version compatibility is not guaranteed, but most normalisation changes etc will cancel. The resulting particle diffusion coefficient are calculated using \cite Lyons1974A \cite Lyons1974B, \cite Albert2005 and such. Note that this makes no sense for E fields!
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions. A set of test arguments is supplied. Call using ./calculate_diffusion `<test_pars` to use these. Or try ./calculate_diffusion -h for argument help
\verbinclude help_i.txt
  \author Heather Ratcliffe \date 17/02/2017
*/

const char PER_UTIL_HELP_ID = 'i';/**<ID to identify help file for this utility*/

struct diff_cmd_line{
  std::string file_prefix;/**< Prefix part of file names*/
  size_t d[2];/**< Dimensions of D to produce*/
  bool is_list;/**< Use FFT or spectrum list input*/
  bool is_spect;/**< Use spectrum list*/
  size_t ref;/**< Reference sdf file number to get Bx info from*/
  std::string ref_name;/**< Reference file name to use for Bx info*/
  std::vector<std::string> file_list;
  int fuzz;/**<Fuzz for spectral cutout*/
  int smth;/**<Smoothing width for output spectrum*/
  size_t n_ang;/**<Number of angles for output spectrum*/
  int wave;/**<Wave type ID (see support.h WAVE_* )*/
  int ang;/**< Angular function type (can be FUNCTION_NULL) */
};/**< \brief Command line arguments for diffusion calculation utility*/

diff_cmd_line special_command_line(int argc, char *argv[]);
bool is_filenumber(char * str);

/** \brief Main program
*
* Calculate particle diffusion from a spectrum
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char *argv[]){

  bool err;
  size_t per_proc = 1;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);
  
  MPI_Barrier(MPI_COMM_WORLD);

  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
//  setup_args cmd_line_args = process_command_line(argc, argv);
  diff_cmd_line cmd_line_args = special_command_line(argc, argv);

  if(cmd_line_args.is_list){
    per_proc = std::ceil( (float) cmd_line_args.file_list.size() / (float) mpi_info.n_procs);
  }

  /* Get constants from deck and share to other procs*/
  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  log_code_constants(cmd_line_args.file_prefix);

  my_print("Processing "+mk_str(per_proc)+" blocks per core", mpi_info.rank);

  data_array Bx;
  if(cmd_line_args.ref_name ==""){
    reader my_reader = reader(cmd_line_args.file_prefix);
    my_reader.update_ref_filenum(cmd_line_args.ref);
    size_t my_space[2];
    my_space[0] = extract_space_part(cmd_line_args.file_list[0]).first;
    my_space[1] = extract_space_part(cmd_line_args.file_list[cmd_line_args.file_list.size()-1]).second;
    Bx = get_Bx(cmd_line_args.file_prefix, my_space, cmd_line_args.ref);
  }else{
    Bx = data_array(cmd_line_args.file_prefix +cmd_line_args.ref_name);
  }

  controller contr = controller(cmd_line_args.file_prefix);

  //---------------- Now we loop over blocks per proc-------
  for(size_t block_num = 0; block_num < per_proc; block_num++){
    data_array dat_fft;
    calc_type spec_norm = 1.0;
    std::string filename;
    if(per_proc*mpi_info.rank + block_num < cmd_line_args.file_list.size()){
      filename = cmd_line_args.file_prefix + cmd_line_args.file_list[per_proc*mpi_info.rank + block_num];

      //Check version compatibility.
      if(!check_wipps_version(filename)){
        //If breaking changes made, add an exit at this point.
      }
      size_t n_tims = extract_num_time_part(filename);
      if(!cmd_line_args.is_spect){
        dat_fft = data_array(filename);
        if(!dat_fft.is_good()){
          my_error_print("Failed to read "+filename);
          continue;
        }
        //Have an FFT, need a spectrum
        contr.set_plasma_B0(dat_fft.B_ref);
        //Match spectrum omega size to FFT omega size
        size_t spec_sz = dat_fft.get_dims(0);
/** \todo generate_spectrum goes wrong if spec_sz != k_sz Note we have to rebin if we change to work*/
        contr.add_spectrum(spec_sz, cmd_line_args.n_ang, (cmd_line_args.ang != FUNCTION_NULL));

        //Renorm. FFT
        if(dat_fft.space[1] > 1) spec_norm = 0.5 * dat_fft.space[1]*dat_fft.space[1];
        if(n_tims > 0) spec_norm *= n_tims;
        //Spectrum is actually B^2 so we square the FFT norming factor
        spec_norm *= spec_norm;
        //2D data and want to use angle function, so squash in k_y
        if(dat_fft.get_dims() == 3 && cmd_line_args.ang != FUNCTION_NULL){
          dat_fft = dat_fft.total(1);
        }
        contr.get_current_spectrum()->generate_spectrum(dat_fft, cmd_line_args.fuzz, cmd_line_args.ang);
      }else{
        //Read spectrum file
        err = contr.add_spectrum(filename);
        if(err){
          my_error_print("Failed to read "+filename);
          continue;
        }
        if(contr.get_current_spectrum()->space[1] > 1) spec_norm = 0.5 * std::pow(contr.get_current_spectrum()->space[1], 2);
        if(n_tims > 0) spec_norm *= n_tims;
      /** \todo We have to do the time bins better than this! */

        spec_norm *= spec_norm;
      }
      calc_type om_ce = contr.get_plasma().get_omega_ref("ce");
    /** \todo Add truncation limits from cmd line */
      if(spec_norm > 0){
        contr.get_current_spectrum()->apply(spectrum::part::B, divide, spec_norm);
      }else{
        my_error_print("Erroneous spectrum normalisation constant. Proceeeding unnormalised");
      }
    }else{
      continue;
      //We've run out of files, must have been non-integer decomposition
    }
    //Now we have some spectral data, from an FFT dump or from a spectrum dump
    contr.add_d(cmd_line_args.d[0], cmd_line_args.d[1]);
    contr.get_current_d()->calculate();

  }
  //-----------------end of per_proc loop---- Now controller holds one spectrum and d per block

  MPI_Barrier(MPI_COMM_WORLD);//Wait until all blocks are complete
  
  bounce_av_data bounce_dat;

  //Set up our bounce-averaging params and do the average
  //Accumulated property is for every file in a directory
  bounce_dat.set_Bx_size(Bx.get_dims(0));
  bounce_dat.set_Bx(Bx);
  bounce_dat.max_latitude = 90.0;
  bounce_dat.L_shell = 4.5;
  contr.bounce_average(bounce_dat);
  
  contr.save_spectra(cmd_line_args.file_prefix);
  contr.save_D(cmd_line_args.file_prefix);

  MPI_Finalize();
  //call these last...

  exit(0);
}

bool is_filenumber(char * str){
//Filenumber would be positive int without any leading spaces

  for(size_t i=0; i< strlen(str); i++){
    if(str[i] < 0 || str[i] > 9) return false;
  }
  return true;
}

diff_cmd_line special_command_line(int argc, char *argv[]){
/** \brief Process commandline arguments
*
* This handles all command line arguments to this utility, so expects no arguments not listed below or in spect_process_command_line
*/

  diff_cmd_line values;
  //Default values
  values.file_prefix = "./files/";
  values.ref = 0;
  values.ref_name = "";
  values.d[0] = 50;
  values.d[1] = 50;
  values.is_list = false;
  values.is_spect = false;
  values.file_list.clear();

  spect_args spec_vals = spect_process_command_line(argc, argv);
  values.fuzz = spec_vals.fuzz;
  values.smth = spec_vals.smth;
  values.n_ang = spec_vals.n_ang;
  values.wave = spec_vals.wave;
  values.ang = spec_vals.ang;

  for(int i = 1; i < argc; i++){
    if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-d")==0 && i < argc-2){
      if(atoi(argv[i+1]) >= 0 && atoi(argv[i+2]) >= 0){
        values.d[0] = atoi(argv[i+1]);
        values.d[1] = atoi(argv[i+2]);
      }else{
        my_error_print("D not an integer or is negative!");
      }
      i+=2;
    }
    else if(strcmp(argv[i], "-ref")==0 && i < argc-1){
      if(is_filenumber(argv[i+1])) values.ref = atoi(argv[i+1]);
      else values.ref_name = argv[i+1];
      i++;
    }
    else if(((strcmp(argv[i], "-Finput")==0)||(strcmp(argv[i], "-Sinput")==0)) && i < argc-1){
      values.is_list=true;
      if((strcmp(argv[i], "-Sinput")==0)) values.is_spect = true;
      //Now hunt for next arg..., we assume no '-' starting filenames. The filenames will be extracted later
      int tmp = i;
      while(i<argc-1 && argv[i+1][0]!= '-' && argv[i+1][0]!=HANDLED_ARG[0]) i++;
      if(tmp -i >= 1 ) i--;
      //Go back one so that loop advance leaves us in correct place, but not if we didn't skip on at all or we'd infinite loop
    }else if(!((strlen(argv[i]) > 0) && argv[i][0] == HANDLED_ARG[0])){
      std::cout<<"UNKNOWN OPTION " <<argv[i]<<'\n';
    }
  }
  if(values.d[0] >MAX_SIZE){
    values.d[0] = MAX_SIZE;
    my_error_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  if(values.d[1] >MAX_SIZE){
    values.d[1] = MAX_SIZE;
    my_error_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }

  if(!values.is_list) my_error_print(std::string("Specify either -Sinput or -Finput"));

    //If using a file-list extract the names
    if(values.is_list){
      values.file_list = process_filelist(argc, argv);
    }
  return values;
}
/** @} */



