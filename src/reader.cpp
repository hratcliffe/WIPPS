//
//  reader.cpp
//  
//
//  Created by Heather Ratcliffe on 02/10/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "support.h"
#include "reader.h"
#include "my_array.h"
#include "sdf.h"
#include <mpi.h>



/** This is going to take care of openeing and reading from SDF. We'll give it a pointer to a storage class instance, and it can then use that to obtain pointer to the bit of memory to fill. We'll also give it a numerical range, a file prefix, and a block name to grab. Perhaps we also have option to restrict on the ranges, so we can easily block in space.
*
*It will be sort of safe. It wont store the data array to write to, but will need to be fed it. So we can;t overwrite the memory.
*
*/

reader::reader(std::string file_prefix_in,  char * block_id_in){
/**We'll also have to work out how many zeros to use at some point. For nw we'll do it when we make the reader eh, as that's where we assign the prefix. We'll assume a 0th file exists, and default and minimum is 4. We'll also try 5-7.
*/
  strcpy(this->block_id, block_id_in);
  this->file_prefix = file_prefix_in;
  
  std::ifstream file;
  std::string name = file_prefix + "000.sdf";
  n_z=3;

  while(!file.is_open()){
   name.insert(file_prefix.size(), "0");
   file.open(name);
   n_z ++;
   if(n_z > 8) break;
  }
  file.close();

}

bool reader::read_dims(int &n_dims, std::vector<int> &dims){

//opens 0th file, gets dimension info so we can construct an array.... Note we don't have to read the data, only the block list.
  char fmt[5];
  char file_num[10];
  sprintf(fmt,"%s%d%c" , "%0", n_z, 'd');
  snprintf(file_num, 10, fmt, 0);
  std::string file_name = file_prefix + file_num +".sdf";

  std::cout<<"Opening "<<file_name<<std::endl;
  
  sdf_file_t *handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block;

  if(handle){std::cout<<"Success! File opened"<<std::endl;}
  else{
    std::cout<<"Bleh. File open error. Aborting"<<std::endl;
    return 1;
  }

  bool err=sdf_read_blocklist(handle);
  if(!err){std::cout<<"Yay! Blocks read"<<std::endl;}
  else{std::cout<<"Bum. Block read failed"<<std::endl;}


  block = sdf_find_block_by_id(handle, this->block_id);

  if(block){std::cout<<"Success!!! Requested block found"<<std::endl;}
  else{
    std::cout<<"Bleh. Requested block not found. Aborting"<<std::endl;
    return 1;
  }

  n_dims = block->ndims;
  for(int i=0; i<n_dims; i++) dims.push_back(block->dims[i]);
  
  sdf_close(handle);

  return 0;
}

bool reader::read_data(data_array * my_data_in, int time_range[2], int space_range[2]){
/** This will open the files dictated by time range seuqentially, and populate them into the data_array. It'll stop when the end of rnage is reached, or it goes beyond the size available. Space range upper entry of -1 is taken as respective limit. @return 0 for success, 1 for error \todo Currently gives no report of nature of error...
*/

//First we check the space range is within range, and chnage the -1's to the repsective dimensions. If out of rnage we return error (1).
  int dim = my_data_in->get_dims(0);
  if((space_range[1] > dim) || ((space_range[1]> 0)  &&(space_range[0] > space_range[1]))) return 1;

  if(space_range[0]==-1) space_range[0] = 0;
  if(space_range[1]==-1) space_range[1] = dim;
  
  std::cout<<space_range[0]<<" "<<space_range[1]<<std::endl;

  if(time_range[0] < 0) time_range[0] = 0;
  if(time_range[1] < time_range[0]) time_range[1] = time_range[0];

  strcpy(my_data_in->block_id, block_id);
  //set block id

  //Now we start from time_range[0] and run through files to time_range[1], or until file not found. We construct with sprintf and the known n_z

  char fmt[5];
  char file_num[10];
  sprintf(fmt,"%s%d%c" , "%0", n_z, 'd');
  std::string file_name;
  sdf_file_t *handle;
  sdf_block_t * block;
  bool err;
  my_type * ax_ptr;
  int len;

  snprintf(file_num, 10, fmt, time_range[0]);
  file_name = file_prefix + file_num +".sdf";

  std::cout<<"Opening "<<file_name<<std::endl;
  
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);

  //first we open first file and do grids.
  block = sdf_find_block_by_id(handle, "grid");
  handle->current_block = block;

  //checks type is plain i/e/ even gridded mesh
  if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){std::cout<< "Uh oh, grids look wrong"<<std::endl;}

  sdf_read_data(handle);

  ax_ptr = my_data_in->get_axis(0, len);
  memcpy ((void *)ax_ptr, block->grids[0], len*sizeof(my_type));
  /**get 0th axis. \todo extend to 2-d data */
  sdf_close(handle);

  ax_ptr = my_data_in->get_axis(my_data_in->get_dims()-1, len);
  //pointer to last axis, which will be time

  int i;
  //now loop over files and get actual data
  for(i=time_range[0]; i<time_range[1];++i){

    snprintf(file_num, 10, fmt, i);
    file_name = file_prefix + file_num +".sdf";

    std::cout<<"Opening "<<file_name<<std::endl;
  
    handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
    if(!handle) break;

    sdf_read_blocklist(handle);
    block = sdf_find_block_by_id(handle, this->block_id);
    if(!block) break;
  
    handle->current_block = block;
    sdf_read_data(handle);

    *(ax_ptr + i) = handle->time;
    //save time of file
    if(!block->data) break;

    my_data_in->populate_row(block->data, dim, i-time_range[0]);

    sdf_close(handle);

  }

//report if we broke out of loop and print filename

  if(i < time_range[1]-1){
    std::cout<<"Read stopped by error at file "<<file_name<<std::endl;
    return 1;
  }

return 0;
}

