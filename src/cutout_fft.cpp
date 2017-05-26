
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

/** \ingroup utils
 */

/** \defgroup cutout_util FFT cutout utility
*@{ *\brief Utility to trim FFTd data to specified axis limits
*
*Reads array from given file, cuts out to supplied limits and saves to given output file. deck.status file is read from file_prefix+deck.status and allows to specify frequency cuts in w_ce. Wavenumber cuts are in m^-1. If no output file is given, output will be in [inputfile]_trim. Output file gains current version number, so we enforce strict checking
\verbinclude help_c.txt
\author Heather Ratcliffe \date 11/08/2016

*/

const char PER_UTIL_HELP_ID = 'c';/**<ID to identify help file for this utility*/

struct cutout_args{
  std::string file_in;/**< Input FFT file*/
  std::string file_out;/**< Output FFT file*/
  std::string file_prefix;/**<File path prepended to all filenames*/
  std::vector<my_type> limits;/**<Limits to use for cutout*/

};/**< \brief Command line arguments for FFT cutout utility*/

cutout_args cutout_process_command_line(int argc, char *argv[]);

/** \brief Main program
*
* Read an FFT, cutout to limits and print to new file
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char *argv[]){

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
  
  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  cutout_args my_args = cutout_process_command_line(argc, argv);

  if(mpi_info.rank == 0){

    get_deck_constants(my_args.file_prefix);
    

    //Check version compatibility.
    if(!check_wipps_version(my_args.file_prefix+my_args.file_in)){
      //Exit at this point, since output file will inherit new version specifier but probably wont be right
      exit(1);
    }
    
    data_array FFT_in = data_array(my_args.file_prefix+my_args.file_in);
    //Set cutout limits on FFT
    size_t n_dims = FFT_in.get_dims();
    
    if(my_args.limits.size() !=2*n_dims){
      my_error_print("******Please supply 2 limits per dimension*****", mpi_info.rank);
      exit(1);
    }
    if(my_const.omega_ce != 0){
      size_t sz = my_args.limits.size();
      
      my_args.limits[sz-2] *= my_const.omega_ce;
      my_args.limits[sz-1] *= my_const.omega_ce;
      
    }
    //If we've got a deck.status, we assume last 2 limits are in terms of omega_ce
    
    std::fstream file;
    file.open((my_args.file_prefix + my_args.file_out).c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()){
      FFT_in.write_section_to_file(file, my_args.limits);
    }
    file.close();
    my_print( "FFT section output in "+(my_args.file_prefix + my_args.file_out), mpi_info.rank);
  }

  MPI_Finalize();
  //call these last...
  
  exit(0);
}

cutout_args cutout_process_command_line(int argc, char *argv[]){
/** \brief Process commandline arguments
*
* Expects full list and no more.
*/

  cutout_args values;
  values.file_prefix = "./files/";
  values.file_in = "";
  values.file_out = "";

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
    }
    else if(strcmp(argv[i], "-lims")==0 && i < argc-1){
      while(i<argc-1 && (argv[i+1][0]!= '-'  || ((argv[i+1][1] >='0' && argv[i+1][1] <='9') || argv[i+1][1] =='.'))){
        //Checks if next argument is a new flag, but allows negative numbers
        values.limits.push_back(checked_strtof(argv[i+1]));
        i++;
      }
    }
    else my_error_print(std::string("UNKNOWN OPTION ")+argv[i]);
    i++;
  }
  if(values.file_out == "" && values.file_in != "") values.file_out = append_into_string(values.file_in, "_trim");
  
  return values;

}
/** @} */
/** @} */

