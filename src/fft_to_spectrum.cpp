//
//  FFT_to_spectrum.cpp
//  
//
//  Created by Heather Ratcliffe on 12/08/2016.
//
//

#include <stdio.h>
#include <cstring>
#include "support.h"
#include "data_array.h"
#include "controller.h"
#include "spectrum.h"

/** \ingroup utils */

/** \defgroup spect_util Spectrum generation utility
*@{ 
*\brief Utility to generate a spectrum from Fourier transformed data
*
*Requires an input directory and an input fft'd data file. The output is either specified, or is the input file with _spectrum appended before the extension. The "wave" option specifies the wave mode by single-character key (w, p, o) and defaults to Whistler. A "fuzz" parameter controlling how tight a band around the dispersion curve can be supplied as a percentage, default is 10%. Spectra contain both frequency and angle data, the n_ang, ang and extra flags control this. FFTs may not stay compatible cross-code version, so we do a version check first.
\verbinclude help_f.txt
\author Heather Ratcliffe \date 12/08/2016
\todo Add per-util target to makefile
*/

const char PER_UTIL_HELP_ID = 'f';/**<ID to identify help file for this utility*/

struct fft_spect_args{
  std::string file_prefix;/**< Filepath prepended to all files*/
  std::string file_in;/**<Input FFT filename*/
  std::string file_out;/**<Output spectrum filename (optional)*/
  int fuzz;/**<Fuzz for spectral cutout*/
  int smth;/**<Smoothing width for output spectrum*/
  size_t n_ang;/**<Number of angles for output spectrum*/
  float ang_sd;/**< Width for angular function (if applicable)*/
  int wave;/**<Wave type ID (see support.h WAVE_* )*/
  int ang;/**< Angular function type (can be FUNCTION_NULL) */
  bool mask;/**< Flag to output spectrum extraction mask to file also*/
};/**< \brief Command line arguments for spectrum generation utility*/


fft_spect_args fft_spect_process_command_line(int argc, char *argv[]);

/** \brief Main program
*
* Convert FFT to spectrum and write to file
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char *argv[]){
/** \todo FFT normalisation -> V2.0. Create file converter if so*/
//We don't need MPI here but SDF does

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);

  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  fft_spect_args my_args = fft_spect_process_command_line(argc, argv);
  
  get_deck_constants(my_args.file_prefix);

  log_code_constants(my_args.file_prefix);
  /* Do this so we have the spectrum wavetype and angle codes*/

  controller contr = controller(my_args.file_prefix);
    //For spectrum we need a plasma to dictate dispersion, so we go via controller

  //Check version compatibility.
  if(!check_wipps_version(my_args.file_prefix+my_args.file_in)){
    //If breaking changes made, add an exit at this point.
  }
  
  data_array data_in = data_array(my_args.file_prefix+my_args.file_in);

  if(!data_in.is_good()){
      my_error_print("Data array allocation failed. Aborting.");
      return 0;
  }

  my_type B_ref = 0.0;
  if(std::abs(data_in.B_ref) > 1e-12) B_ref = data_in.B_ref;
  else B_ref = my_const.omega_ce * me/std::abs(q0);
  //Use reference B from file if non-zero, else use deck constants

  size_t k_len = data_in.get_dims(0);//K_x length
  contr.add_spectrum(k_len, my_args.n_ang, (my_args.ang != FUNCTION_NULL));

  if(data_in.get_dims() == 3 && my_args.ang != FUNCTION_NULL){
    //2D data and want to use angle function, so squash in k_y
    data_in = data_in.total(1);
  }
  contr.set_plasma_B0(B_ref);

  data_array * mask;
  if(my_args.mask){
    //Mask requested, create array for it and call
    mask = new data_array();
    mask->clone_empty(data_in);
    contr.get_current_spectrum()->generate_spectrum(data_in, my_args.fuzz, my_args.ang, my_args.ang_sd, mask);
  }else{
    contr.get_current_spectrum()->generate_spectrum(data_in, my_args.fuzz, my_args.ang, my_args.ang_sd);
  }
  //Do smoothing
  if(my_args.smth > 1) contr.get_current_spectrum()->smooth_B(my_args.smth);
  //Write output
  std::fstream outfile;
  outfile.open(my_args.file_prefix+my_args.file_out, std::ios::binary|std::ios::out);
  contr.get_current_spectrum()->write_to_file(outfile);
  outfile.close();

  //Write mask output if required
  if(my_args.mask){
    std::string filename =my_args.file_prefix+my_args.file_out;
    filename = append_into_string(filename, "_mask");
    outfile.open(filename, std::ios::binary|std::ios::out);
    mask->write_to_file(outfile);
    outfile.close();

  }
  
  return 0;
}


fft_spect_args fft_spect_process_command_line(int argc, char *argv[]){
/** \brief Process special command line args
*
*Process the fft utility arguments. Expects full list and no more.
\todo Change to a protected integer conversion
*/

  fft_spect_args values;
  //Default values if nothing supplied
  values.file_prefix = "./files/";
  values.file_in = "";
  values.file_out = "";
 
  spect_args spec_vals = spect_process_command_line(argc, argv);
  values.fuzz = spec_vals.fuzz;
  values.smth = spec_vals.smth;
  values.mask = spec_vals.mask;
  values.n_ang = spec_vals.n_ang;
  values.wave = spec_vals.wave;
  values.ang = spec_vals.ang;
  values.ang_sd = spec_vals.ang_sd;
  
  int i = 1;
  while(i < argc){
    if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-in")==0 && i < argc-1){
      values.file_in = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-out")==0 && i < argc-1){
      values.file_out = argv[i+1];
      i++;
    }else if(!((strlen(argv[i]) > 0) && argv[i][0] == HANDLED_ARG[0])){
      std::cout<<"UNKNOWN OPTION " <<argv[i]<<'\n';
    }
    i++;
  }
  if(values.file_out == "" && values.file_in != "") values.file_out = append_into_string(values.file_in, "_spectrum");
  //catch accidental overwriting request
  if(values.file_out == values.file_in) values.file_out = append_into_string(values.file_in, "_spectrum");
  return values;

}
/** @} */
/** @} */
