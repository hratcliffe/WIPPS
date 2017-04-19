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

/** \defgroup utils Utility programs
*@{ */

/** \defgroup spect_util Spectrum generation utility
*@{ 
*\brief Utility to generate a spectrum from Fourier transformed data
*
*Requires an input directory and an input fft'd data file. The output is either specified, or is the input file with _spectrum appended before the extension. The "wave" option specifies the wave mode by single-character key (w, p, o) and defaults to Whistler. A "fuzz" parameter controlling how tight a band around the dispersion curve can be supplied as a percentage, default is 10%. Spectra contain both frequency and angle data, the n_ang, ang and extra flags control this. 
\verbinclude help_f.txt
\author Heather Ratcliffe \date 12/08/2016
*/

const char PER_UTIL_HELP_ID = 'f';/**<ID to identify help file for this utility*/

struct fft_spect_args{
  std::string file_prefix;/**< Filepath prepended to all files*/
  std::string file_in;/**<Input FFT filename*/
  std::string file_out;/**<Output spectrum filename (optional)*/
  int fuzz;/**<Fuzz for spectral cutout*/
  int smth;/**<Smoothing width for output spectrum*/
  size_t n_ang;/**<Number of angles for output spectrum*/
  int wave;/**<Wave type ID (see support.h WAVE_* )*/
  int ang;/**< Angular function type (can be FUNCTION_NULL) */
  bool mask;/**< Flag to output spectrum extraction mask to file also*/
};/**< Command line arguments for this utility*/


fft_spect_args fft_spect_process_command_line(int argc, char *argv[]);

int main(int argc, char *argv[]){
/** \todo FFT normalisation -> V2.0*/
//We don't need MPI here but SDF does

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);

  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  fft_spect_args my_args = fft_spect_process_command_line(argc, argv);
  
  get_deck_constants(my_args.file_prefix);

  log_code_constants(my_args.file_prefix);
  /* Do this so we have the spectrum wavetype and angle codes*/

  controller contr = controller(my_args.file_prefix);
    //For spectrum we need a plasma to dictate dispersion, so we go via controller
  data_array data_in = data_array(my_args.file_prefix+my_args.file_in, 0);

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
    contr.get_current_spectrum()->generate_spectrum(data_in, my_args.fuzz, my_args.ang, DEFAULT_SPECTRUM_ANG_STDDEV, mask);
  }else{
    contr.get_current_spectrum()->generate_spectrum(data_in, my_args.fuzz, my_args.ang);
  }
  //Do smoothing
  if(my_args.smth >1) contr.get_current_spectrum()->smooth_B(my_args.smth);
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
*/

  fft_spect_args values;
  //Default values if nothing supplied
  values.file_prefix = "./files/";
  values.file_in = "";
  values.file_out = "";
  values.fuzz = 10;
  values.smth = 0;
  values.mask = false;
  values.n_ang = DEFAULT_N_ANG;
  values.wave = WAVE_WHISTLER;
  values.ang = FUNCTION_DELTA;
  
  bool extr = false;
  
  for(int i=1; i< argc; i++){
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
    }
    else if(strcmp(argv[i], "-om")==0 && i < argc-1){
      values.fuzz = atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-n_ang")==0 && i < argc-1){
      values.n_ang = atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-wave")==0 && i < argc-1){
      if(argv[i+1][0] == 'w' || argv[i+1][0] == 'W') values.wave=WAVE_WHISTLER;
      else if(argv[i+1][0] == 'p' || argv[i+1][0] == 'P') values.wave=WAVE_PLASMA;
      else if(argv[i+1][0] == 'o' || argv[i+1][0] == 'O') values.wave=WAVE_O;
      i++;
    }
    else if(strcmp(argv[i], "-ang")==0 && i < argc-1 && !extr){
      if(argv[i+1][0] == 'd' || argv[i+1][0] == 'D') values.ang=FUNCTION_DELTA;
      else if(argv[i+1][0] == 'g' || argv[i+1][0] == 'G') values.ang=FUNCTION_GAUSS;
      else if(argv[i+1][0] == 'i' || argv[i+1][0] == 'I') values.ang=FUNCTION_ISO;
      i++;
    }
    else if(strcmp(argv[i], "-extr")==0){
      values.ang = FUNCTION_NULL;
      extr = true;
      //This _overrides_ -ang
    }
    else if(strcmp(argv[i], "-mask")==0){
      values.mask = true;
    }
    else if(strcmp(argv[i], "-smooth")==0 && i < argc-1){
      values.smth = atoi(argv[i+1]);
      i++;
    }else{
      std::cout<<"UNKNOWN OPTION " <<argv[i]<<'\n';
    }
  }
  if(values.file_out == "" && values.file_in != "") values.file_out = append_into_string(values.file_in, "_spectrum");
  //catch accidental overwriting request
  if(values.file_out == values.file_in) values.file_out = append_into_string(values.file_in, "_spectrum");
  if(values.smth < 0) values.smth = 0;
  return values;

}
/** @} */
/** @} */
