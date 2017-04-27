//
//  compress_distributions.cpp
//  

#include <stdio.h>
#include <math.h>
#include "main.h"
#include "support.h"
#include "reader.h"
#include "data_array.h"

/** \ingroup utils */
/** \defgroup dist_util Distribution extraction utility
*@{ 
*\brief Utility to extract distribution functions from data
*
*Extracts distributions from SDF files and optionally compresses their spatial dimension. 
\verbinclude help_d.txt
\author Heather Ratcliffe \date 09/09/2016
\caveat This utility only works for 2-D distributions right now.
\todo Fix to work for 1, 2, or 3-d? Remember all distribs are 3-d with absent dims of size 1
*/

const char PER_UTIL_HELP_ID = 'd';/**<ID to identify help file for this utility*/

struct dist_cmd_line{
  bool list;/**<Flag to just list the available distribs */
  std::string file_prefix;/**<File path prepended to all files*/
  int dump;/**<Dump number to read*/
  int blocks;/**<Number of space blocks to divide x dimension into*/
};/**< \brief Command line arguments for distribution utility*/


dist_cmd_line special_command_line(int argc, char *argv[]);

/** \brief Main program
*
* Extract distributions and space average
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
  
  process_command_line_help_arg(argc, argv, PER_UTIL_HELP_ID);
  dist_cmd_line my_args = special_command_line(argc, argv);

  char id[ID_SIZE] = "Distrib";
  reader my_reader = reader(my_args.file_prefix, id, my_args.dump);

  std::string special_tag = "dist_fn";
  std::vector<std::pair<std::string, std::string> > dist_blocks;

  {  //Read block list and create list of just the distrib fns
    std::vector<std::pair<std::string, std::string> > blocks = my_reader.list_blocks();
    for(size_t i=0; i< blocks.size(); i++){
      if(special_tag.compare(blocks[i].first.substr(0, special_tag.size())) ==0) dist_blocks.push_back(blocks[i]);
    }
  }
  if(dist_blocks.size() == 0) my_print("No distribution functions found");
  if(my_args.list){
    for(size_t i=0; i< dist_blocks.size(); i++) my_print(dist_blocks[i].first);
    exit(0);
  }
  //Now we extract any distributions and flatten them on x into the specified number of blocks, and output
  data_array dat;
  size_t n_dims;
  std::vector<size_t> dims;
  size_t dims_arr[2];
  for(size_t i=0; i< dist_blocks.size(); i++){
    my_reader.read_dims(n_dims, dims, dist_blocks[i].second);
    for(size_t j=0; j<2; j++) dims_arr[j] = dims[j];
    dat = data_array((size_t)2, dims_arr);
    if(!dat.is_good()){
      my_error_print("Error allocating array for "+dist_blocks[i].second+", skipping");
      continue;
    }
    
    //The distribs are always 3-d, with absent dims of size 1. Ours should be 2-d initially
    err = my_reader.read_distrib(dat, dist_blocks[i].second, my_args.dump);
    if(err){
      my_error_print("Error reading distrib "+dist_blocks[i].second);
      continue;
    }
    std::fstream file;
    size_t trim_size = 8; //Len of "dist_fn/" leading part
    std::string tag = replace_char((dist_blocks[i].first).substr(trim_size, std::string::npos), '/', '_');
    std::string filename = my_args.file_prefix+tag+'_'+mk_str(my_args.dump) +".dat";
    my_print("Writing "+dist_blocks[i].second);
    file.open(filename.c_str(),std::ios::out|std::ios::binary);

    if(my_args.blocks ==-1){
      //Average all blocks
      dat=dat.average(0);
      dat.write_to_file(file, false);
      dat.write_closer(file);
    }
    else if(my_args.blocks ==0){
      //Simply read and then dump
      dat.write_to_file(file, false);
      dat.write_closer(file);
    }
    else{
      data_array section;
      size_t block_size = ceil(dat.get_dims(0) / (float)my_args.blocks);
      for(int j = 0; j< my_args.blocks; j++){
        //Just current block in x
        section=dat.average(0, block_size*j, block_size*(j+1));
        section.space[0] =block_size*j;
        section.space[1] =block_size*(j+1);
        //Close file on final block
        section.write_to_file(file, false);
      }
      //Close the file
      dat.write_closer(file);
    }

    my_print("Closing");
    file.close();
  }

}

dist_cmd_line special_command_line(int argc, char *argv[]){
/** \brief Process commandline arguments
*
* This handles all command line arguments to this utility, so expects no arguments not listed below
*/


  dist_cmd_line values;
  //Default values
  values.list=0;
  values.file_prefix = "./files/";
  values.dump = 0;
  values.blocks = -1;

  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-list")==0){
      values.list = 1;
    }
    else if(strcmp(argv[i], "-extr")==0){
      values.blocks = 0;
    }
    else if(strcmp(argv[i], "-dump")==0 && i < argc-1){
      values.dump= atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-x_blocks")==0 && i < argc-1){
      values.blocks= atoi(argv[i+1]);
      i++;
    }
    else my_error_print(std::string("UNKNOWN OPTION ")+argv[i]);
    
  }
  return values;
}
/** @} */
/** @} */

