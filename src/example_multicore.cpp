
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

/** \ingroup example_mains */

const char PER_UTIL_HELP_ID = 'm';/**<ID to identify help file for this utility*/

struct my_args{
  std::string file_prefix;/**< Prefix part of file names (path plus any common prefix) */
  int num_in;/**< An integer*/
  std::vector<std::string> file_list;/** A list of files to work on*/
};/**< \brief Command line arguments for example utility*/

my_args example_process_command_line(int argc, char *argv[]);

/** \brief Main program
*
* Do things
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char *argv[]){

  //Initialise MPI
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }
  
  //Pass the rank: by default this function prints only from root (rank = 0)
  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  //This checks whether the '-h' flag was given, and if so prints the help text from file help_[x] where x=PER_UTIL_HELP_ID
  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);

  //This does our other command line stuff, filling a my_args struct with the results
  my_args cmd_line_args = example_process_command_line(argc, argv);

  //Calculate how many files per processor to work on
  size_t per_proc = std::ceil( (float) cmd_line_args.file_list.size() / (float) mpi_info.n_procs);
  MPI_Barrier(MPI_COMM_WORLD);
  //Print from all cores, because of the -1 argument here
  my_print("This is processor " + mk_str(mpi_info.rank), mpi_info.rank, -1);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);

  std::string file_in = "";
   //--------------Loop over blocks per proc-------
  for(size_t file_num = 0; file_num < per_proc; file_num++){

    file_in = cmd_line_args.file_list[file_num];
    //Check version compatibility of input data file. If anything is changed that affects how code should operate, this can be used to adapt for it
    if(!check_wipps_version(cmd_line_args.file_prefix + file_in)){
      //exit(1);//Uncomment to exit at this point if required
      //continue;//Uncomment to go to next file now if required
    }
    
    //Read from file_in into a data array
    data_array input_data = data_array(cmd_line_args.file_prefix + file_in);
    
    size_t n_dims = input_data.get_dims();
    //Make a copy of the data
    data_array output_data = input_data;
    //////Do some stuff to the array
    if(output_data.is_good()){
      //Divide elementwise by 2
      output_data.apply(divide, 2);
      //Subtract  input_data from output_data elementwise
      output_data.apply(subtract, input_data);
      //Average on 0th dim, if one exists
      if(n_dims > 0) output_data.average(0);

      //////Write data to file
      std::fstream file;
      //Create a filename by putting "_out" into the input name
      std::string filename;
      filename = cmd_line_args.file_prefix + file_in;
      //Add the rank to filename to ensure we don't accidentally write to the same file
      filename = append_into_string(filename, "_out" + mk_str(mpi_info.rank));
      //Open file for output
      file.open(filename.c_str(),std::ios::out|std::ios::binary);
      if(file.is_open()){
        output_data.write_to_file(file);
        file.close();
        my_print( "Data output in "+filename, mpi_info.rank, -1);
      }
    }
  }

  //Clean-up MPI
  MPI_Finalize();
  //call these last...
  
  exit(0);
}

my_args example_process_command_line(int argc, char *argv[]){
/** \brief Process commandline arguments
*
* Expects full list and no extra options, so will warn if it doesn't understand a key
*/

  my_args values;
  values.file_prefix = "./files/";
  values.num_in = 0;
  values.file_list.clear();
  
  int i = 1;
  while(i < argc){
    if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }else if(strcmp(argv[i], "-num")==0 && i < argc-1){
      values.num_in = checked_strtol(argv[i+1]);

    }else if(strcmp(argv[i], "-files")==0 && i < argc-1){
      while(i<argc-1 && argv[i+1][0]!= '-'){
        //Checks if next argument is a new flag
        values.file_list.push_back(argv[i+1]);
        i++;
      }
    }
    else my_error_print(std::string("UNKNOWN OPTION ")+argv[i], mpi_info.rank);
    i++;
  }
  return values;

}
/** @} */

