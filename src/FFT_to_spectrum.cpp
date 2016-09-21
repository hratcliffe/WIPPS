//
//  FFT_to_spectrum.cpp
//  
//
//  Created by Heather Ratcliffe on 12/08/2016.
//
//
/** \file FFT_to_spectrum.cpp \brief Utility to convert an FFT to a spectrum
*
* Reads supplied FFT file and cuts out spectrum. Assumes Whistler mode
*/


#include <stdio.h>

#include "support.h"
#include "main_support.h"
#include "data_array.h"
#include "controller.h"
#include "spectrum.h"

extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/
deck_constants my_const;/**< Physical constants*/

struct FFT_spect_args{
  std::string file_prefix;
  std::string file_in;
  std::string file_out;
  int fuzz;
  size_t n_ang;
  int wave;
  int ang;
};

FFT_spect_args FFT_spect_process_command_line(int argc, char *argv[]);


int main(int argc, char *argv[]){
/** \todo FFT normalisation*/
//We don't need MPI here but SDF does

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  
  FFT_spect_args my_args = FFT_spect_process_command_line(argc, argv);
  
  get_deck_constants(my_args.file_prefix);

  log_code_constants(my_args.file_prefix);
  /* Do this so we have the spectrum wavetype and angle codes*/

  controller contr = controller(my_args.file_prefix);
  //For spectrum we need a plasma to dictate dispersion, so we go via controller
  
  data_array data_in = data_array(my_args.file_prefix+my_args.file_in, 0, false);

  if(!data_in.is_good()){
    my_print("Data array allocation failed. Aborting.");
    return 0;
  }

  my_type B_ref = 0.0;
  if(std::abs(data_in.B_ref) > 1e-12) B_ref = data_in.B_ref;
  else B_ref = my_const.omega_ce * me/std::abs(q0);
  //Use reference B from file if non-zero, else use deck constants

  size_t k_len = data_in.get_dims(0);//K_x length
  contr.add_spectrum(k_len, my_args.n_ang, (my_args.ang != FUNCTION_NULL));

  contr.set_plasma_B0(B_ref);
  contr.get_current_spectrum()->generate_spectrum(data_in, my_args.fuzz, my_args.ang);

  data_array BB = contr.get_current_spectrum()->copy_out_B();
  std::cout<<BB.maxval()<<'\n';
  
  std::fstream outfile;
  outfile.open(my_args.file_prefix+my_args.file_out, std::ios::binary|std::ios::out);
  contr.get_current_spectrum()->write_to_file(outfile);
  outfile.close();
  
  
}

FFT_spect_args FFT_spect_process_command_line(int argc, char *argv[]){

  FFT_spect_args values;
  values.file_prefix = "./files/";
  values.file_in = "";
  values.file_out = "";
  values.fuzz = 10;
  values.n_ang = DEFAULT_N_ANG;
  values.wave = WAVE_WHISTLER;
  values.ang = FUNCTION_DELTA;
  
  bool extr = false;
  
  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-h")==0){
      print_help('f');
      exit(0);
    }
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
    }else{
      std::cout<<"UNKNOWN OPTION " <<argv[i]<<'\n';
    }
  }
  if(values.file_out == "" && values.file_in != "") values.file_out = append_into_string(values.file_in, "_spectrum");
  
  return values;

}
