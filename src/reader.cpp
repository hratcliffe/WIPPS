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
#include <algorithm>
#include "support.h"
#include "main_support.h"
#include "reader.h"
#include "data_array.h"
#include "sdf.h"
#include <mpi.h>


extern const mpi_info_struc mpi_info;

reader::reader(){
  n_z =0;
  ref_file_num = 0;
  file_prefix="";

  time_range[0]=0; time_range[1]=0; time_range[2]=0;
  space_range[0]=0; space_range[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  
}
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

bool reader::read_dims(size_t &n_dims, std::vector<size_t> &dims){
/** \brief Gets dimensions of the block specified in reader.
*
*Opens reference file, and gets dimension info. Returns by reference, with 0 for success, 1 for file open or read failure. Note we don't have to read the data, only the block list.
*/

  std::string file_name = get_full_name(ref_file_num);

  my_print("Getting dimensions", mpi_info.rank);
  sdf_file_t *handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block;

  if(!handle){
    my_print("Bleh. File open error. Aborting", mpi_info.rank);
    return 1;
  }
  
  bool err=sdf_read_blocklist(handle);
  if(err) my_print("Block read failed", mpi_info.rank);

  block = sdf_find_block_by_id(handle, this->block_id);

  if(!block){
    my_print("Requested block not found. Aborting", mpi_info.rank);
    return 1;
  }

  if(block->datatype != my_sdf_type) my_print("WARNING!!! Data type does not match. Output may be nonsense!", mpi_info.rank);

  if(!is_accum(this->block_id)){
    n_dims = block->ndims;
  }else{
    n_dims = block->ndims -1 ;
    //remove time dimension
  }
  for(size_t i=0; i<n_dims; i++) dims.push_back(block->dims[i]);
  
  sdf_close(handle);

  return 0;
}

int reader::read_data(data_array &my_data_in, int time_range[3], int space_range[2], int flatten_on){
/** \brief Read data into given array
*
*This will open the files dictated by time range sequentially, and populate them into the data_array. It'll stop when the end of range is reached, or it goes beyond the size available. Space range upper entry of -1 is taken as respective limit. NB: blocking is only supported on the X axis. @return 0 for success, 1 for error 2 for unusual exit, i.e. early termination \todo Break out into three functions, common, plain and acc
*/
  
  if(!my_data_in.is_good()){
    my_print("Cannot read into invalid array", mpi_info.rank);
    return 1;
  }
  strcpy(my_data_in.block_id, block_id);
  //set block id

  //Now we start from time_range[0] and run through files to time_range[1], or until file not found.

  sdf_file_t *handle;
  sdf_block_t * block ,*ax_block;
  my_type * ax_ptr;
  size_t len;

  std::string file_name = get_full_name(time_range[0]);
  bool accumulated = is_accum(block_id);

  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);

  std::string grid_name = "grid";
  if(accumulated) grid_name = "grid_accum";
    
  //first we open first file and do grids.
  block = sdf_find_block_by_id(handle, grid_name.c_str());
  handle->current_block = block;

  //checks type is plain i/e/ even gridded mesh
  if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){
    my_print("Uh oh, grids look wrong", mpi_info.rank);
  }
  sdf_read_data(handle);

  if(block->datatype != my_sdf_type){
    //Axis type doesn't match expected, we'll get nonsense
    my_print("Wrong data type detected. Converting grids", mpi_info.rank);
  }

  for(size_t i=0, i2=0; i< my_data_in.get_dims()-1; i++, i2++){
    if(flatten_on >= 0 && flatten_on == i){
      i2++;
      continue;
    }
    //Skip a flattening dim if exists
    ax_ptr = my_data_in.get_axis(i, len);

    if(block->datatype != my_sdf_type) std::copy((other_type *) block->grids[i2], (other_type *) block->grids[i2] + len, ax_ptr);
    else std::copy((my_type *) block->grids[i2], (my_type *) block->grids[i2] + len, ax_ptr);
    //We assume if the type does not match then it's the other of double or float
    //Get space axes
  }
  bool simple_slice=false;
  size_t * source_sizes = nullptr;
  int source_advance = 1;
  //Simple slices are those where we're slicing only the last spatial dimension. For what we have here, that means 1 space dim or more than one and no slicing. We have no mechanism to slice y etc, only x.
  //Grid dimensions are for staggered grid for subtract 1
  if(my_data_in.get_dims() ==2 || block->dims[0]-1 == space_range[1]-space_range[0]){
  //Either result is 2-d, i.e. single time data is 1-d, or we're not trying to crop 0th dim
    simple_slice = true;
  }
  source_sizes = (size_t *) malloc((block->ndims)*sizeof(size_t));
  
  for(size_t i=0; i< block->ndims - accumulated; i++){
    //Get SPACE advance for source, i.e advance between one time read and the next
    source_sizes[i] = block->dims[i]-1;
    //Minus one because this is grid block with stagger
    source_advance*=source_sizes[i];

  }
  my_type * flat_data = nullptr;
  if(flatten_on >=0){
    //One time at once
    flat_data = (my_type*) malloc(sizeof(my_type)*source_advance);
  }

  sdf_close(handle);

  ax_ptr = my_data_in.get_axis(my_data_in.get_dims()-1, len);
  //pointer to last axis, which will be time
  int i;
  int last_report=0;
  int report_interval = (time_range[1]-time_range[0])/10;
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;
    //Say we want to report 10 times over the list, or every 20th file if  more than 200.

  my_data_in.space[0] = space_range[0];
  my_data_in.space[1] = space_range[1];
  
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
    
    my_ptr += source_advance/block->dims[0]*space_range[0];

    if(!accumulated){
      size_t val[1];
      val[0] = i-time_range[0];
      *(ax_ptr + val[0]) = (my_type) handle->time;

      if(simple_slice && flatten_on < 0){
        //Simple slicing and no flattening required
        my_data_in.populate_slice(my_ptr, 1, val);
      }else if(simple_slice){
        //Simple slicing but we flatten first
        flatten_fortran_slice(my_ptr, flat_data, my_data_in.get_dims()-1, source_sizes,flatten_on);
        my_data_in.populate_slice(flat_data, 1, val);
      }else{
        my_data_in.populate_complex_slice(my_ptr, 1, val, source_sizes);
      }
      total_reads++;
    }
    else{
      ax_block = sdf_find_block_by_id(handle, grid_name.c_str());
      handle->current_block = ax_block;
      sdf_read_data(handle);

      //don't read more than time[2] rows
      rows = ax_block->dims[block->ndims-1];
      if(total_reads + rows >= time_range[2]) rows = time_range[2]- total_reads;

      //Copy time grid out
      if(ax_ptr){
        if(ax_block->datatype != my_sdf_type) std::copy((other_type *) ax_block->grids[1], (other_type *) ax_block->grids[1] + rows, ax_ptr+total_reads);
        else std::copy((my_type *) ax_block->grids[1], (my_type *) ax_block->grids[1] + rows, ax_ptr+total_reads);
      }
      size_t val[1];
      if(simple_slice && flatten_on < 0){
        for(int j=0; j<rows; j++){
          val[0] = total_reads+j;
          my_data_in.populate_slice(my_ptr, 1, val);
          my_ptr += source_advance;
        }
      }else if(simple_slice){
        for(int j=0; j<rows; j++){
          val[0] = total_reads+j;
          //Simple slicing but we flatten first
          flatten_fortran_slice(my_ptr, flat_data, my_data_in.get_dims()-1, source_sizes,flatten_on);

          my_data_in.populate_slice(flat_data, 1, val);
          my_ptr += source_advance;
        }
        
      }else{
        for(int j=0; j<rows; j++){
          val[0] = total_reads+j;
          my_data_in.populate_complex_slice(my_ptr, 1, val, source_sizes);
          my_ptr += source_advance;
        }
      }
      total_reads+= rows;
    }
    sdf_close(handle);
    if(accumulated && total_reads >=time_range[2]) break;
  }

  if(source_sizes) free(source_sizes);

  //Set 0th time
  size_t n_dims = my_data_in.get_dims();
  my_data_in.time[0] = my_data_in.get_axis_element(n_dims-1, 0);
  
  //report if we broke out of loop and print filename
  if(i < time_range[1] || (accumulated && total_reads < time_range[2])){
    my_data_in.resize(1, total_reads);
    //trim array to number of lines read

    //Set last time
    my_data_in.time[1] = my_data_in.get_axis_element(n_dims-1, my_data_in.get_dims(n_dims-1)-1);
    
    if(!accumulated) my_print("Read stopped by error at file "+file_name, mpi_info.rank);
    else my_print("Read "+mk_str(total_reads)+" rows", mpi_info.rank);
    return 2;
  }
  else{
    //Set last time
    my_data_in.time[1] = my_data_in.get_axis_element(n_dims-1, my_data_in.get_dims(n_dims-1)-1);
    return 0;
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
  int sz =block->next_block_location;
  sdf_close(handle);

  return sz;


}

bool reader::current_block_is_accum(){
  return is_accum(this->block_id);
}
bool reader::is_accum(std::string block_id){
/** \brief checks for time accumulated blocks */

  if(block_id == "ax" || block_id =="ay" || block_id =="az" || block_id =="abx" || block_id =="aby" || block_id =="abz") return true;


  return false;

}
