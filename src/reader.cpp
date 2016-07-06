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


extern const mpi_info_struc mpi_info;

reader::reader(std::string file_prefix_in,  char * block_id_in, int first){
/** \brief Create reader
*
*Sets up ids and tests how many zeros are in filename by trial and error. First is by default 0, and is the reference dump number to use for testing dims and zeros.
Default and minimum is 4 or the length of first as a string, tries 5-7 also.
*/
  strcpy(this->block_id, block_id_in);
  this->file_prefix = file_prefix_in;
  this->ref_file_num = first;
  std::ifstream file;

  std::string file_ref = mk_str(ref_file_num);
  int n_z_first = file_ref.size();
  if(n_z_first > MAX_FILENAME_DIGITS){
    my_print("Really that many files??", 0, mpi_info.rank);
  }

  std::string name = file_prefix + file_ref+".sdf";
  n_z=n_z_first;

  while(!file.is_open()){
    name.insert(file_prefix.size(), "0");
    file.open(name);
    n_z ++;
    if(n_z > MAX_FILENAME_DIGITS+1) break;
  }
  file.close();

}

bool reader::read_dims(int &n_dims, std::vector<int> &dims){
/** \brief Gets dimensions of the block specified in reader.
*
*Opens reference file, and gets dimension info. Returns by reference, with 0 for success, 1 for file open or read failure. Note we don't have to read the data, only the block list.

*/

  std::string file_name = get_full_name(ref_file_num);

  my_print("Getting dimensions", mpi_info.rank);
  sdf_file_t *handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block;

  if(handle){my_print("Success! File opened", mpi_info.rank);}
  else{
    my_print("Bleh. File open error. Aborting", mpi_info.rank);
    return 1;
  }

  bool err=sdf_read_blocklist(handle);
  if(!err) my_print("Yay! Blocks read", mpi_info.rank);
  else my_print("Bum. Block read failed", mpi_info.rank);

  block = sdf_find_block_by_id(handle, this->block_id);

  if(block) my_print("Success!!! Requested block found", mpi_info.rank);
  else{
    my_print("Bleh. Requested block not found. Aborting", mpi_info.rank);
    return 1;
  }

  if(block->datatype != my_sdf_type) my_print("WARNING!!! Data type does not match. Output may be nonsense!", mpi_info.rank);

  if(!is_accum(this->block_id)){
    n_dims = block->ndims;
  }else{
    n_dims = block->ndims -1 ;
    //remove time dimension
  }
  for(int i=0; i<n_dims; i++) dims.push_back(block->dims[i]);
  
  sdf_close(handle);

  return 0;
}

int reader::read_data(data_array * my_data_in, int time_range[3], int space_range[2]){
/** \brief Read data into given array
*
*This will open the files dictated by time range sequentially, and populate them into the data_array. It'll stop when the end of range is reached, or it goes beyond the size available. Space range upper entry of -1 is taken as respective limit. @return 0 for success, 1 for error \todo Currently gives no report of nature of error... use 2 for non-fatal read error \todo 2-d array???!
*/
  
  strcpy(my_data_in->block_id, block_id);
  //set block id

  //Now we start from time_range[0] and run through files to time_range[1], or until file not found. We construct with sprintf and the known n_z

  sdf_file_t *handle;
  sdf_block_t * block;
  my_type * ax_ptr;
  int len;

  std::string file_name = get_full_name(time_range[0]);
  
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);

  //first we open first file and do grids.
  block = sdf_find_block_by_id(handle, "grid");
  handle->current_block = block;

  //checks type is plain i/e/ even gridded mesh
  if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){
    my_print("Uh oh, grids look wrong", mpi_info.rank);
  }
  sdf_read_data(handle);

  if(block->datatype != my_sdf_type){
    //Axis type doesn't match expected, we'll get nonsense
    my_print("Wrong data type detected. Grids will be corrupt", mpi_info.rank);
  }

  for(int i=0; i< my_data_in->get_dims()-1; i++){
    ax_ptr = my_data_in->get_axis(i, len);

    // Mostly c++ way
    std::copy((my_type *) block->grids[i], (my_type *) block->grids[i] + len, ax_ptr);
    /**get 0th axis. \todo extend to 2-d data */
  }
  sdf_close(handle);

  ax_ptr = my_data_in->get_axis(my_data_in->get_dims()-1, len);
  //pointer to last axis, which will be time
  bool accumulated = is_accum(block_id);
  int i;
  int last_report=0;
  int report_interval = (time_range[1]-time_range[0])/10;
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;
    //Say we want to report 10 times over the list, or every 20th file if  more than 200.

  my_data_in->time[0] = time_range[0];
  my_data_in->time[1] = time_range[1];
  my_data_in->space[0] = space_range[0];
  my_data_in->space[1] = space_range[1];

  int rows=0;
  int total_reads=0;
  //now loop over files and get actual data

  for(i=time_range[0]; i<time_range[1];++i){
    file_name = get_full_name(i);

    if((i-last_report) >= report_interval){
      my_print("Opening " + file_name, mpi_info.rank);

      last_report = i;
    }

    handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
    if(!handle) break;

    sdf_read_blocklist(handle);
    block = sdf_find_block_by_id(handle, this->block_id);
    if(!block) break;
  
    handle->current_block = block;
    sdf_read_data(handle);

    if(block->datatype != my_sdf_type) my_print("WARNING!!! Data type does not match. Output may be nonsense!", mpi_info.rank);

    //save time of file
    if(!block->data) break;
    my_type * my_ptr = (my_type *) block->data;
    my_ptr +=space_range[0];
    
    if(!accumulated){
      *(ax_ptr + i) = (my_type) handle->time;

      my_data_in->populate_row(my_ptr, space_range[1], i-time_range[0]);
      total_reads++;

    }
    else{
      block = sdf_find_block_by_id(handle, "grid_accum");
      handle->current_block = block;
      sdf_read_data(handle);

      rows = block->dims[block->ndims-1];
      if(total_reads + rows >= time_range[2]) rows = time_range[2]- total_reads;
      //don't read more than time[2] rows
      std::cout<<total_reads<<" "<<rows<<" "<<len<<'\n';
      for(int j=0; j< rows; j++) std::cout<< *((my_type*) block->grids[1] + j)<<" ";
      //if(ax_ptr) std::copy((my_type *) block->grids[1], (my_type *) block->grids[1] + rows, ax_ptr +total_reads);
      if(ax_ptr) for(int j=0; j< rows-1; j++) *(ax_ptr+total_reads + j) = *((my_type*) block->grids[1] + j);
      if(ax_ptr) for(int j=0; j< rows; j++) std::cout<< *(ax_ptr+total_reads + j)- *((my_type*) block->grids[1] + j) <<" ";

      //Copy time grid out
      for(int j=0; j<rows; j++){
        my_data_in->populate_row(my_ptr, space_range[1], total_reads+j);
        my_ptr += block->dims[0];
      }
      total_reads+= rows;
    }
    sdf_close(handle);
    if(accumulated && total_reads >=time_range[2]) break;
  }

  //report if we broke out of loop and print filename
  if(i < time_range[1]){
    my_data_in->resize(1, total_reads);
    //trim array to number of lines read... NB 2-D ONLY
    if(!accumulated) my_print("Read stopped by error at file "+file_name, mpi_info.rank);
    else my_print("Read "+mk_str(total_reads)+" rows", mpi_info.rank);
    return 2;
  }

return 0;
}

std::string reader::get_full_name(int num){
/** \brief Construct file name
*
*Combines the prefix, correct number of zeros and .sdf extension into name
*/

  char fmt[5];
  char file_num[10];
  sprintf(fmt,"%s%d%c" , "%0", n_z, 'd');

  snprintf(file_num, 10, fmt, num);
  return file_prefix + file_num +".sdf";

}

int reader::get_file_size(){
/**\brief Get recorded file size
*
*Acts as basic check of SDF integrity and reading. This should match size on disk.
*/
  sdf_file_t *handle;
  sdf_block_t * block;

  std::string file_name = get_full_name(0);
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);

  block = handle->last_block_in_file;
  return block->next_block_location;


}

bool reader::is_accum(std::string block_id){
/** \brief checks for time accumulated blocks */

  if(block_id == "ax" || block_id =="ay" || block_id =="az" || block_id =="abx" || block_id =="aby" || block_id =="abz") return true;


  return false;

}

