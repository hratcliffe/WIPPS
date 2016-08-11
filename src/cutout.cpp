/** \file cutout.cpp \brief Helper program to make a cutout of an array
*
*
* Reads array from given file, cuts out to supplied limits and saves to given output file. deck.status file is read from file_prefix+deck.status and allows to specify frequency cuts in w_ce and space in ... If no output file is given, output will be in [inputfile]_trim
Call example:

\author Heather Ratcliffe \date 11/08/2016.
*/

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <mpi.h>

#include "main.h"
#include "support.h"
#include "reader.h"
#include "data_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "main_support.h"

deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/

struct cutout_args{

  std::string file_in;
  std::string file_out;
  std::string file_prefix;
  std::vector<my_type> limits;

};

cutout_args cutout_process_command_line(int argc, char *argv[]);

int main(int argc, char *argv[]){

  int err;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);
  if(mpi_info.n_procs > 1){
    my_print("Utility uses one core only!", mpi_info.rank);
  }
  
  cutout_args my_args = cutout_process_command_line(argc, argv);

  if(mpi_info.rank == 0){

    get_deck_constants(my_args.file_prefix);
    
    data_array FFT_in = data_array(my_args.file_prefix+my_args.file_in);
    //Set cutout limits on FFT
    int n_dims = FFT_in.get_dims();

    for(int i=0; i< my_args.limits.size(); i++){
      std::cout<<my_args.limits[i]<<' ';
    }
    std::cout<<my_args.file_prefix<<" "<<my_args.file_in<<" "<<my_args.file_out<<'\n';
    
    if(my_args.limits.size() !=2*n_dims){
      my_print("******Please supply 2 limits per dimension*****", mpi_info.rank);
      exit(1);
    }
    if(my_const.omega_ce != 0){
      size_t sz = my_args.limits.size();
      
      my_args.limits[sz-2] *= my_const.omega_ce;
      my_args.limits[sz-1] *= my_const.omega_ce;
      
    }
    //If we've got a deck.status, we assume last 2 limits are in terms of omega_ce
    
    //Construct filename. Since the MPI is using block-wise domain decomposition, different processors can't overlap on blocks
    std::fstream file;
    file.open((my_args.file_prefix + my_args.file_out).c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()){
      FFT_in.write_section_to_file(file, my_args.limits);
    }
    file.close();
    my_print( "FFT section output in "+(my_args.file_prefix + my_args.file_out), mpi_info.rank);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...


  exit(0);
}

cutout_args cutout_process_command_line(int argc, char *argv[]){


  cutout_args values;
  values.file_prefix = "./files/";
  values.file_in = "";
  values.file_out = "";

  for(int i=0; i< argc; i++){
    if(strcmp(argv[i], "-h")==0) print_help('c');
    
    if(strcmp(argv[i], "-f")==0 && i < argc-1) values.file_prefix = argv[i+1];
    if(strcmp(argv[i], "-in")==0 && i < argc-1) values.file_in = argv[i+1];
    if(strcmp(argv[i], "-out")==0 && i < argc-1) values.file_out = argv[i+1];
    if(strcmp(argv[i], "-lims")==0 && i < argc-1){
      while(i<argc-1 && (argv[i+1][0]!= '-'  || ((argv[i+1][1] >='0' && argv[i+1][1] <='9') || argv[i+1][1] =='.'))){
        //Checks if next argument is a new flag, but allows negative numbers
        values.limits.push_back(atof(argv[i+1]));
        i++;
      }
    }
  }
  if(values.file_out == "" && values.file_in != "") values.file_out = append_into_string(values.file_in, "_trim");
  
  return values;

}
