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
*Sets up ids, sets n_z etc. NOTE n_z must be correctly set before any reads are done. Either by supplying reference filenumber (any sdf file which exists) or using update_ref_filenum. 
*/
  strcpy(this->block_id, block_id_in);
  this->file_prefix = file_prefix_in;
  this->ref_file_num = first;
  std::ifstream file;

  std::string file_ref = mk_str(ref_file_num);
  int n_z_first = file_ref.size();
  if(n_z_first > MAX_FILENAME_DIGITS){
    this->ref_file_num = std::pow(10, MAX_FILENAME_DIGITS);
    file_ref=mk_str(ref_file_num);
    n_z_first = file_ref.size();
    my_error_print("Reference file number out of range, truncating to "+file_ref, 0, mpi_info.rank);
  }
  this->n_z = get_filename_n_z(ref_file_num);
}

int reader::get_filename_n_z(int file_num){
/** Check on disk to see how many zeros in filenames by attempting to read a file of number filenum at various lengths. Default and minimum is 4 or the length of first as a string, tries up to MAX_FILENAME_DIGITS.*/

  std::ifstream file;

  std::string file_ref = mk_str(file_num);
  int n_z_first = file_ref.size();

  std::string name = file_prefix + file_ref+".sdf";
  int n_z=n_z_first;
  file.open(name);
  while(!file.is_open()){
    name.insert(file_prefix.size(), "0");
    file.open(name);
    n_z ++;
    if(n_z > MAX_FILENAME_DIGITS+1) break;
  }
  file.close();
  return n_z;

}

std::vector<std::pair<std::string, std::string> > reader::list_blocks(){
/** \brief List blocks in ref file
*
*Returns a vector containing pairs of block name and the internal id.
*/
  std::string file_name = get_full_name(this->ref_file_num);
  sdf_file_t *handle;
  sdf_block_t * next;
  std::vector<std::pair<std::string, std::string> > list;
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle){
    my_error_print("Error reading file "+file_name);
    return list;
  }
  sdf_read_blocklist(handle);

  next = handle->current_block;
  for(int i =0; i< handle->nblocks; i++){
    list.push_back(std::make_pair(std::string(next->name), std::string(next->id)));
    next = next->next;
  }
  sdf_close(handle);
  return list;
}


bool reader::read_dims(size_t &n_dims, std::vector<size_t> &dims){
  return read_dims(n_dims, dims, std::string(this->block_id));
}

bool reader::read_dims(size_t &n_dims, std::vector<size_t> &dims, std::string b_id){
/** \brief Gets dimensions of the block specified in reader.
*
*Opens reference file, and gets dimension info. Returns by reference, with 0 for success, 1 for file open or read failure. Note we don't have to read the data, only the block list.
*/

  std::string file_name = get_full_name(ref_file_num);

  my_print("Getting dimensions", mpi_info.rank);
  sdf_file_t *handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block;

  if(!handle){
    my_error_print("File open error "+file_name, mpi_info.rank);
    return 1;
  }
  
  bool err=sdf_read_blocklist(handle);
  if(err) my_error_print("Block read failed", mpi_info.rank);

  block = sdf_find_block_by_id(handle, b_id.c_str());

  if(!block){
    my_error_print("Requested block not found", mpi_info.rank);
    return 1;
  }

  if(!is_accum(b_id)){
    n_dims = block->ndims;
  }else{
    n_dims = block->ndims -1 ;
    //remove time dimension from accum. blocks
  }
  for(size_t i=0; i<n_dims; i++) dims.push_back(block->dims[i]);
  
  sdf_close(handle);

  return 0;
}

int reader::pre_read(data_array& my_data_in, int ref_time, bool accumulated, int flatten_on, size_t &n_dims, size_t * &source_sizes){
/** \brief Read grids
*
*Reads grids blocks into supplied array and also returns the source sizes for reference.
*/
  sdf_file_t *handle;
  sdf_block_t * block;
  my_type * ax_ptr;
  size_t len;

  std::string file_name = get_full_name(ref_time);
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);

  std::string grid_name = "grid";
  if(accumulated) grid_name = "grid_accum";
    
  //Locate correct grids block
  block = sdf_find_block_by_id(handle, grid_name.c_str());
  handle->current_block = block;

  //checks type is plain i/e/ even gridded mesh
  if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){
    my_error_print("Uh oh, grids look wrong", mpi_info.rank);
  }
  sdf_read_data(handle);

  if(block->datatype != my_sdf_type){
    my_error_print("Wrong data type detected. Converting grids", mpi_info.rank);//Print a warning about needing to type-convert
  }

  //Copy (all) the spatial axes
  for(size_t i=0, i2=0; i< my_data_in.get_dims()-1; i++, i2++){
    if(flatten_on >= 0 && (size_t)flatten_on == i){
      i2++;
      continue;
    }
    //Skip a flattening dim if exists
    ax_ptr = my_data_in.get_axis(i, len);

    if(block->datatype != my_sdf_type) std::copy((other_type *) block->grids[i2], (other_type *) block->grids[i2] + len, ax_ptr);
    else std::copy((my_type *) block->grids[i2], (my_type *) block->grids[i2] + len, ax_ptr);
    //We assume if the type does not match then it's the other of double or float
  }

  //Now we get the *spatial* input sizes
  n_dims = block->ndims - accumulated;
  source_sizes = (size_t *) malloc(n_dims*sizeof(size_t));

  for(size_t i=0; i< n_dims; i++){
    //Size of data source. Note -1 because we're doing this from the grid block
    source_sizes[i] = block->dims[i]-1;
  }

  sdf_close(handle);

  return 0;
}


int reader::read_data(data_array &my_data_in, size_t time_range[3], size_t space_range[2], int flatten_on){
/** \brief Read data into given array
*
*This will open the files dictated by time range sequentially, and populate them into the data_array. It'll stop when the end of range is reached, or it goes beyond the size available. Space range upper entry of -1 is taken as respective limit. NB: blocking is only supported on the X axis. @return 0 for success, 1 for error 2 for unusual exit, i.e. early termination 
*/
  
  if(!my_data_in.is_good()){
    my_error_print("Cannot read into invalid array", mpi_info.rank);
    return 1;
  }
  strcpy(my_data_in.block_id, block_id);
  //set block id

  sdf_file_t *handle;
  sdf_block_t * block;
  std::string file_name;

  bool accumulated = is_accum(block_id);
  if(flatten_on >= (int) my_data_in.get_dims()){
    my_print("Flatten dimension exceeds array size. No flattening will occur", mpi_info.rank);
    flatten_on = -1;
  }
  
  size_t grid_n_dims;
  size_t * source_sizes = nullptr;
  int err = pre_read(my_data_in, time_range[0], accumulated, flatten_on, grid_n_dims, source_sizes);
  //Do the pre-read loops stuff, like grids and size getting
  if(err) return 1;

  bool simple_slice=false;
  int source_advance = 1;
  
  //Simple slices are those where we're slicing only the last spatial dimension, so we can just stop the read early. Since we only slice in x, that means either one and only one spatial dim, or no slicing.
  simple_slice = (my_data_in.get_dims() ==2 || source_sizes[0] == space_range[1]-space_range[0]);

  //Get SPACE advance for source, i.e advance between one time read and the next
  for(size_t i=0; i< grid_n_dims; i++){
    source_advance*=source_sizes[i];
  }

  //If we have to flatten, we need somewhere to put the flat data
  my_type * flat_data = nullptr;
  if(flatten_on >=0){
    flat_data = (my_type*) malloc(sizeof(my_type)*source_advance);
  }

  //pointer to last axis, i.e. the time axis time
  size_t i;
  size_t last_report=0;
  size_t report_interval = (time_range[1]-time_range[0])/10;
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;
    //Say we want to report 10 times over the list, or every 20th file if  more than 200.

  my_data_in.space[0] = space_range[0];
  my_data_in.space[1] = space_range[1];
  
  size_t rows=0;
  size_t total_reads=0;
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

    if(!block->data) break;
    my_type * my_ptr = (my_type *) block->data;
    
    //Move on to the correct 0th dim starting location
    my_ptr += source_advance/block->dims[0]*space_range[0];

    if(!accumulated){
      size_t val[1];
      val[0] = i-time_range[0];
      read_plain_time(my_data_in, handle, val[0]);
      
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
      //Read all rows in file, but not more than time[2] rows total
      rows = block->dims[block->ndims-1];
      if(total_reads + rows >= time_range[2]) rows = time_range[2]- total_reads;

      read_acc_time(my_data_in, handle, total_reads, rows);
      size_t val[1] = {total_reads};
      for(size_t j=0; j<rows; j++){

        if(simple_slice && flatten_on < 0){
          my_data_in.populate_slice(my_ptr, 1, val);
        }else if(simple_slice){
          //Simple slicing but we flatten first
          flatten_fortran_slice(my_ptr, flat_data, my_data_in.get_dims()-1, source_sizes,flatten_on);
          my_data_in.populate_slice(flat_data, 1, val);
        }else{
          my_data_in.populate_complex_slice(my_ptr, 1, val, source_sizes);
        }
        val[0] ++;
        my_ptr += source_advance;
      }
      total_reads+= rows;
    }
    sdf_close(handle);
    if(accumulated && total_reads >=time_range[2]) break;
  }

  //If we did flattening, what we actually want is the average not the total
  if(flatten_on >= 0){
    my_data_in.divide(source_sizes[flatten_on]);
  }

  if(source_sizes) free(source_sizes);
  if(flat_data) free(flat_data);

  //Set 0th time
  size_t n_dims = my_data_in.get_dims();
  my_data_in.time[0] = my_data_in.get_axis_element(n_dims-1, 0);

  //report if we broke out of loop and print info
  if((!accumulated && i < time_range[1]) || (accumulated && total_reads < time_range[2])){
    my_data_in.resize(1, total_reads);
    //trim array to number of lines read

    //Set last time
    my_data_in.time[1] = my_data_in.get_axis_element(n_dims-1, my_data_in.get_dims(n_dims-1)-1);
    
    if(!accumulated) my_error_print("Read stopped by error at file "+file_name + " ("+mk_str(total_reads)+" times)", mpi_info.rank);
    else my_print("Read "+mk_str(total_reads)+" times", mpi_info.rank);
    return 2;
  }
  else{
    //Set last time
    my_data_in.time[1] = my_data_in.get_axis_element(n_dims-1, my_data_in.get_dims(n_dims-1)-1);
    my_print("Read "+mk_str(total_reads)+" times", mpi_info.rank);
    return 0;
  }

return 0;
}

int reader:: read_plain_time(data_array& my_data_in, sdf_file_t *handle, size_t pos){
/** Read time for non-accumulated data*/
  size_t len;
  my_type * ax_ptr = my_data_in.get_axis(my_data_in.get_dims()-1, len);
  *(ax_ptr + pos) = (my_type) handle->time;
  
  return 0;
}

int reader::read_acc_time(data_array & my_data_in, sdf_file_t * handle, size_t total_reads, size_t rows){
/** \brief Read accumulator time axis into my_data_in
*
*
*/
  sdf_block_t * ax_block = sdf_find_block_by_id(handle, "grid_accum");
  handle->current_block = ax_block;
  sdf_read_data(handle);

  size_t len;
  my_type * ax_ptr = my_data_in.get_axis(my_data_in.get_dims()-1, len);

  size_t n_grids = ax_block->ndims;
  //Copy time grid out
  if(ax_ptr){
    if(ax_block->datatype != my_sdf_type) std::copy((other_type *) ax_block->grids[n_grids-1], (other_type *) ax_block->grids[n_grids-1] + rows, ax_ptr+total_reads);
    else std::copy((my_type *) ax_block->grids[n_grids-1], (my_type *) ax_block->grids[n_grids-1] + rows, ax_ptr+total_reads);
  }
  return 0;
}


bool reader::read_distrib(data_array & my_data_in, std::string dist_id, int dump_number){
/** \brief Read distribution function
*
*Reads the distribution function dist_id into data_in. ID has dist_fn removed and is trimmed to 10 chars max. If data type does not match compiled type we cast to match. 
*/
  if(!my_data_in.is_good()){
    my_error_print("Cannot read into invalid array", mpi_info.rank);
    return 1;
  }
  sdf_file_t *handle;
  sdf_block_t * block;
  std::string file_name = get_full_name(dump_number);

  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);
  block = sdf_find_block_by_id(handle, dist_id.c_str());
  if(!block) return 1;

  handle->current_block = block;
  sdf_read_data(handle);

  if(!block->data) return 1;

  size_t total_els = my_data_in.get_total_elements();
  if(block->datatype != my_sdf_type) my_data_in.populate_data((my_type*)block->data, total_els, 1);
  else my_data_in.populate_data( (my_type*)block->data, total_els, 0);

  strncpy(my_data_in.block_id, dist_id.c_str(), ID_SIZE-1);
  my_data_in.time[0] = handle->time;
  my_data_in.time[1] = handle->time;
  
  
  //Finally read the axes
  block = sdf_find_block_by_id(handle, ("grid/"+dist_id).c_str());
  if(!block) return 1;

  handle->current_block = block;
  sdf_read_data(handle);

  my_type * ax_ptr;
  size_t len;
  for(size_t i=0; i< my_data_in.get_dims(); i++){
    ax_ptr = my_data_in.get_axis(i, len);

    if(block->datatype != my_sdf_type) std::copy((other_type *) block->grids[i], (other_type *) block->grids[i] + len, ax_ptr);
    else std::copy((my_type *) block->grids[i], (my_type *) block->grids[i] + len, ax_ptr);
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
