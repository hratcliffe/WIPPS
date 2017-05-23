//
//  extract_field.cpp
//  
//
//  Created by Heather Ratcliffe on 23/05/2017.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <mpi.h>

#include "support.h"
#include "reader.h"
#include "data_array.h"


/** \ingroup utils */

/** \defgroup extr_util Field extractor utility
*@{ 
*\brief Utility to read sdf files and extract single fields (at first time if accumulated)
*
*Opens files, extracts specified field. Can optionally flatten the data in space.
\verbinclude help_e.txt
  \author Heather Ratcliffe \date 23/05/2017

*/

const char PER_UTIL_HELP_ID = 'e';/**<ID to identify help file for this utility*/

struct extractor_args{
  std::string file_out;/**< Output file*/
  int flat_dim;/**<Dimension to average on*/
};/**< \brief Command line arguments for field extractor utility*/

extractor_args extractor_process_command_line(int argc, char *argv[]);


/** \brief Main program
*
* Extract distributions and space average
  @param argc Command line argument count
  @param argv Command line arguments
  @return System error code */

int main(int argc, char ** argv){

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
  extractor_args extra_args = extractor_process_command_line(argc, argv);
  setup_args cmd_line_args = process_command_line(argc, argv);

  if(mpi_info.rank == 0){

    get_deck_constants(cmd_line_args.file_prefix);
    reader my_reader = reader(cmd_line_args.file_prefix);

    size_t my_space[2];

    size_t times[3] = {cmd_line_args.time[0], cmd_line_args.time[0] + 1, 1};
    //use specified file and read one row, ignoring any additional stuff
    my_reader.update_ref_filenum(times[0]);
    my_reader.change_block_id(cmd_line_args.block);

    size_t n_dims;
    std::vector<size_t> dims;
    int err = my_reader.read_dims(n_dims, dims);
    if(!err && n_dims > 0){
      my_space[0] = 0;
      my_space[1] = dims[0];

      //Create array of right rank and size. Note 0-len dims are not counted in rank
      dims.push_back(1);
      while(dims.size() < 4){dims.push_back(0);}
      data_array data = data_array(dims[0], dims[1], dims[2], dims[3]);
      my_reader.read_data(data, times, my_space);

      //Flatten data if requested
      if(extra_args.flat_dim >= 0 && data.get_dims() > 1){
        data = data.average(extra_args.flat_dim);
      }

      std::fstream file;
      std::string filename;
      if(extra_args.file_out == ""){
        filename = "field_"+cmd_line_args.block+"_"+mk_str(times[0])+".dat";
      }else filename = extra_args.file_out;
      filename = cmd_line_args.file_prefix + filename;

      file.open(filename.c_str(),std::ios::out|std::ios::binary);
      if(file.is_open()){
        data.write_to_file(file);
      }
      file.close();
      my_print("Field "+std::string(cmd_line_args.block)+ " output in "+filename, mpi_info.rank);
    }else{
      my_error_print("Failed to read file");
    }
  }
}

extractor_args extractor_process_command_line(int argc, char *argv[]){
/** \brief Process special arguments to extract_field
*
*Part of argument handling is shared with calculate_diffusion and generate_ffts so we handle only the extras here and must pass the rest on. So we nullify those we handle here to enable warning for unknown arguments by setting the first character of them all to HANDLED_ARG
*/

  extractor_args values;
  //Default values
  values.flat_dim = -1;
  values.file_out = "";
  
  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-flatten")==0 && i < argc-1){
      values.flat_dim = atoi(argv[i+1]);
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      i++;
    }else if(strcmp(argv[i], "-out")==0 && i < argc-1){
      values.file_out = argv[i+1];
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      i++;
    }
  }
  
  return values;
}
/** @} */

