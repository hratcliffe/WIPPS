//
//  compress_distributions.cpp
//  
//
//  Created by Heather Ratcliffe on 09/09/2016.
//
//
//Tiny halper function to total the x_px distributions in x either overall or into n blocks

#include <stdio.h>
#include "main.h"
#include "support.h"
#include "reader.h"
#include "data_array.h"
#include "main_support.h"


struct dist_cmd_line{
  bool list;
  std::string file_prefix;
  int dump;
  
};
deck_constants my_const;/**< Physical constants*/

dist_cmd_line special_command_line(int argc, char *argv[]);

int main(int argc, char *argv[]){

  int err;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }
  dist_cmd_line my_args = special_command_line(argc, argv);

  char id[ID_SIZE] = "Distrib";
  reader my_reader = reader(my_args.file_prefix, id, my_args.dump);

  std::string special_tag = "dist_fn";
  std::vector<std::pair<std::string, std::string> > dist_blocks;

  {  //Read block list and create list of just the distrib fns
    std::vector<std::pair<std::string, std::string> > blocks = my_reader.list_blocks();
    for(int i=0; i< blocks.size(); i++){
      if(special_tag.compare(blocks[i].first.substr(0, special_tag.size())) ==0) dist_blocks.push_back(blocks[i]);
    }
  }
  if(my_args.list){
    for(int i=0; i< dist_blocks.size(); i++) std::cout<<dist_blocks[i].first<<'\n';
    exit(0);
  }
  
  //Now we extract any distributions and flatten them on x into the specified number of blocks, and output
  data_array dat;
  size_t n_dims;
  std::vector<size_t> dims;
  size_t dims_arr[2];
//  for(int i=0; i< dist_blocks.size(); i++){
{  int i=0;
    my_reader.read_dims(n_dims, dims, dist_blocks[i].second);
    for(int j=0; j<2; j++) dims_arr[j] = dims[j];
    dat = data_array((size_t)2, dims_arr);
    /** \todo Why is there a 3rd dim of len 1?*/
    err = my_reader.read_distrib(dat, dist_blocks[i].second, my_args.dump);
    if(err){
      std::cout<<"Error reading distribs"<<'\n';
    }
  
    dat=dat.average(0);
    std::fstream file;
    std::string filename = my_args.file_prefix+"dist_test.dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary);

    dat.write_to_file(file);
    file.close();
  }
  
  
  
}

dist_cmd_line special_command_line(int argc, char *argv[]){
//Do special command line processing here

  dist_cmd_line values;
  values.list=0;
  values.file_prefix = "./files/";
  values.dump = 0;

  for(int i=0; i< argc; i++){
    if(strcmp(argv[i], "-h")==0){
      print_help('d');
      exit(0);
    }
    else if(strcmp(argv[i], "-f")==0 && i < argc-1) values.file_prefix = argv[i+1];
    else if(strcmp(argv[i], "-list")==0){
      values.list = 1;
    }
    else if(strcmp(argv[i], "-dump")==0 && i < argc-1){
      values.dump= atoi(argv[i+1]);
    }

  }
  return values;
}
