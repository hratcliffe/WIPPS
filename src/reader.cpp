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
#include "reader.h"
#include "data_array.h"
#include "sdf.h"
#include <mpi.h>


extern const mpi_info_struc mpi_info;

/********Basic setup and allocation functions ****/

reader::reader(){
/** \brief Create empty reader
*
* Set fields to default null values
*/
  n_z = 0;
  ref_file_num = 0;
  file_prefix = "";

  time_range[0]=0; time_range[1]=0; time_range[2]=0;
  space_range[0]=0; space_range[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  
}

reader::reader(std::string file_prefix_in,  char * block_id_in, int first){
/** \brief Create reader
*
*Sets up ids, sets n_z etc. NOTE n_z must be correctly set before any reads are done. Either by supplying reference filenumber (any sdf file which exists) here or using update_ref_filenum later.
*/
  //Set up some generic things
  time_range[0]=0; time_range[1]=0; time_range[2]=0;
  space_range[0]=0; space_range[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));

  //Setup the specifics
  strncpy(this->block_id, block_id_in, ID_SIZE);
  this->file_prefix = file_prefix_in;
  
  //Make sure reference filenum is in valid range and setup n_z
  first = check_file_num(first);
  this->ref_file_num = first;
  this->n_z = get_filename_n_z(ref_file_num);
}

/********Filename and file manipulators ****/

int reader::check_file_num(int file_num){
/** \brief Get file_num into correct range
*
*Check if file_num is representable in number of chars available. Return it if so, else return 0
*/
  //Stringify and count digits
  std::string file_ref = mk_str(file_num);
  int n_z_first = file_ref.size();
  if(n_z_first > MAX_FILENAME_DIGITS){
    my_error_print("Reference file number out of range, using 0", 0, mpi_info.rank);
    file_num = 0;
  }
  return file_num;
}

int reader::get_filename_n_z(int file_num){
/** \brief Infer number of padding 0s in filenames
*
*Check on disk to see how many zeros in filenames by attempting to read a file of number filenum at various lengths. Default and minimum is 4 or the length of file_num as a string, tries up to MAX_FILENAME_DIGITS.*/

  file_num = check_file_num(file_num);
  std::string file_ref = mk_str(file_num);
  //Pad to minimum size of 4
  if(file_ref.size() < 4) file_ref = std::string(4-file_ref.size(), '0') + file_ref;
  int n_z = file_ref.size();

  //Attempt to read files, inserting '0' padding until one can be opened or we reach MAX_FILENAME_DIGITS
  std::string name = file_prefix + file_ref+".sdf";
  std::ifstream file;
  file.open(name);
  while(!file.is_open()){
    name.insert(file_prefix.size(), "0");
    file.open(name);
    n_z ++;
    if(n_z > MAX_FILENAME_DIGITS) break;
  }
  //Max one file can have been opened so we close it
  file.close();
  return n_z;
}

std::string reader::get_full_name(int file_num){
/** \brief Construct file name
*
*Combines the prefix, correct number of zeros and .sdf extension into name
*/

  //Create proper format string to zero pad filenumber
  //We know n_z is at most MAX_FILENAME_DIGITS, so 2 chars is enough for it
  char fmt[5];
  sprintf(fmt,"%s%d%c" , "%0", n_z, 'd');

  //Create the number string
  char file_num_str[MAX_FILENAME_DIGITS+1];
  snprintf(file_num_str, MAX_FILENAME_DIGITS, fmt, file_num);
  return file_prefix + file_num_str +".sdf";
}

void reader::update_ref_filenum(int num){
/** \brief Update reference file number
*
* Set new reference file number and update n_z parameter*/

  num = check_file_num(num);
  if(num >= 0) this->ref_file_num = num;
  this->n_z = get_filename_n_z(num);
}

int reader::get_file_size(){
/**\brief Get recorded file size
*
*Reads the final block_end position from file. Acts as basic check of SDF integrity and reading. This should match size on disk.
*/
  sdf_file_t * handle;
  sdf_block_t * block;

  std::string file_name = get_full_name(0);
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;

  sdf_read_blocklist(handle);
  block = handle->last_block_in_file;
  int sz = block->next_block_location;
  sdf_close(handle);

  return sz;
}

/********Block manipulators ****/
std::vector<std::pair<std::string, std::string> > reader::list_blocks(){
/** \brief List blocks in reference file
*
*Uses file given by ref_file_num and returns a vector containing pairs of block name and the internal id string.
*/

  //Open file given by ref_file_num
  std::string file_name = get_full_name(this->ref_file_num);
  sdf_file_t * handle;
  sdf_block_t * next;
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  
  //List of blocks in file
  std::vector<std::pair<std::string, std::string> > list;

  if(!handle){
    my_error_print("Error reading file "+file_name);
    return list;
  }
  sdf_read_blocklist(handle);
  
  //Step through blocks, recording the name and id of each
  next = handle->current_block;
  for(int i =0; i< handle->nblocks; i++){
    list.push_back(std::make_pair(std::string(next->name), std::string(next->id)));
    next = next->next;
  }
  sdf_close(handle);
  return list;
}

bool reader::current_block_is_accum(){
/** \brief Check if current block is accumulated
*
*Checks whether the block named in reader is accumulated
*/
  return is_accum(this->block_id);
}

bool reader::is_accum(std::string block_id){
/** \brief Checks for time accumulated blocks
*
* These are an extension of mine to EPOCH proper and have following names. a* are E fields, ab* are B fields.
*/

  if(block_id == "ax" || block_id =="ay" || block_id =="az" || block_id =="abx" || block_id =="aby" || block_id =="abz") return true;
  return false;

}

/********File read routines ****/

int reader::pre_read(data_array& data_in, int ref_time, bool accumulated, int flatten_on, size_t &n_dims, size_t * &source_sizes){
/** \brief Read grids and source sizes
*
* For filenumber ref_time, reads grids blocks into supplied array and also returns the source sizes for reference. The data can be flattened on one dimension before populating array, which can be useful for large data sets as only one copy of the full data will be in use.
  @param my_data_in Data array to fill @param ref_time Filenumber to read @param accumulated Whether block to read is accumulated @param flatten_on Dimension to flatten data on upon reading @n_dims Returns the rank of array read @param source_sizes Returns the dimensions of the data source
*/

  sdf_file_t * handle;
  sdf_block_t * block;
  my_type * ax_ptr;
  size_t len;

  std::string file_name = get_full_name(ref_time);
  handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);
  if(!handle) return 1;
  sdf_read_blocklist(handle);

  //Locate correct grids block by name
  std::string grid_name = "grid";
  if(accumulated) grid_name = "grid_accum";
  block = sdf_find_block_by_id(handle, grid_name.c_str());
  handle->current_block = block;

  //Checks type is plain i.e. an even-gridded mesh
  if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){
    my_error_print("Uh oh, grids look wrong", mpi_info.rank);
  }
  sdf_read_data(handle);

  if(block->datatype != my_sdf_type){
    my_error_print("Wrong data type detected. Converting grids", mpi_info.rank);//Print a warning about needing to type-convert
  }

  //Copy (all) the spatial axes, skipping a flattening dim if one exists
  for(size_t i=0, i2=0; i< data_in.get_dims()-1; i++, i2++){
    if(flatten_on >= 0 && (size_t)flatten_on == i){
      i2++;
      continue;
    }
    ax_ptr = data_in.get_axis(i, len);

    //We assume if the type does not match then it's the other of double or float and convert in the copy
    if(block->datatype != my_sdf_type) std::copy((other_type *) block->grids[i2], (other_type *) block->grids[i2] + len, ax_ptr);
    else std::copy((my_type *) block->grids[i2], (my_type *) block->grids[i2] + len, ax_ptr);
  }

  //Now we get the *spatial* input sizes
  n_dims = block->ndims - accumulated;
  source_sizes = (size_t *) malloc(n_dims*sizeof(size_t));

  for(size_t i=0; i< n_dims; i++){
    //Size of data source. Note -1 because we're doing this from the grid block sizes which have ghost cells
    source_sizes[i] = block->dims[i]-1;
  }

  sdf_close(handle);

  return 0;
}

int reader:: read_plain_time(data_array& data_in, sdf_file_t *handle, size_t pos){
/** \brief Read time for non-accumulated data
*
*Reads time-field from the given sdf filehandle and stores into position pos of data_in time axis
*/

  //Time is always last dimension
  data_in.set_axis_element(data_in.get_dims()-1, pos, (my_type) handle->time);
  return 0;
}

int reader::read_acc_time(data_array & data_in, sdf_file_t * handle, size_t start_pos, size_t n_times){
/** \brief Read time for accumulated data
*
*Reads time axis from the given sdf filehandle and stores n_times values into time axis of my_data starting at start_pos. Converts the data type if needed, assuming it is one of float or double.
*/
  sdf_block_t * ax_block = sdf_find_block_by_id(handle, "grid_accum");
  handle->current_block = ax_block;
  sdf_read_data(handle);

  size_t len;
  my_type * ax_ptr = data_in.get_axis(data_in.get_dims()-1, len);

  size_t n_grids = ax_block->ndims;
  //Copy time grid from ax_block into dat_in.axes
  if(ax_ptr){
    if(ax_block->datatype != my_sdf_type) std::copy((other_type *) ax_block->grids[n_grids-1], (other_type *) ax_block->grids[n_grids-1] + n_times, ax_ptr+start_pos);
    else std::copy((my_type *) ax_block->grids[n_grids-1], (my_type *) ax_block->grids[n_grids-1] + n_times, ax_ptr+start_pos);
  }
  return 0;
}

bool reader::read_dims(size_t &n_dims, std::vector<size_t> &dims){
/** \brief Read dimensions of current block (this->block_id). See reader::read_dims(size_t &n_dims, std::vector<size_t> &dims, std::string b_id)
*/
  return read_dims(n_dims, dims, std::string(this->block_id));
}

bool reader::read_dims(size_t &n_dims, std::vector<size_t> &dims, std::string b_id){
/** \brief Gets dimensions of the block specified by b_id
*
*Looks up block b_id in file numbered ref_file_num and gets dimension info. Returns values in the dims vector parameter and error code 0 for success, 1 for a file open or read failure. Note we don't have to read the data, only the block list.
*/

  //Open file
  std::string file_name = get_full_name(ref_file_num);
  my_print("Getting dimensions", mpi_info.rank);
  sdf_file_t *handle = sdf_open(file_name.c_str(), MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block;

  if(!handle){
    my_error_print("File open error "+file_name, mpi_info.rank);
    return 1;
  }
  
  //Attempt to locate block
  bool err=sdf_read_blocklist(handle);
  if(err){
    my_error_print("Block read failed", mpi_info.rank);
    return 1;
  }
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
  //Copy the dimensions, omitting time if block is accumulated
  for(size_t i=0; i<n_dims; i++) dims.push_back(block->dims[i]);
  
  sdf_close(handle);

  return 0;
}


int reader::read_data(data_array &my_data_in, size_t time_range[3], size_t space_range[2], int flatten_on){
/** \brief Read data into given array
*
*Open files dictated by time_range sequentially, and populates the data_array. Data_array should be set to correct dimensions already, else we return with error and leave data_array in partially updated state. Time_range[0,1] are the file numbers to read, [3] is a number of times (rows) for accumulated data and is ignored for normal blocks. Reading stops when time_range[2] is reached (plain blocks), file number time_range[2] or row time_range[3] is reached (accumulated blocks), or no more files are available on disk. Space range can be set to cut out a part of the x-dimension, or set to the entire x_size. A flattening dimension, flatten_on can be given, in which case @return 0 for success, 1 for error 2 for unusual exit, i.e. early termination
*/
  
  if(!my_data_in.is_good()){
    my_error_print("Cannot read into invalid array", mpi_info.rank);
    return 1;
  }
  
  //Set block id
  strcpy(my_data_in.block_id, block_id);

  sdf_file_t * handle;
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
  
  //Check spatial dims match sliced size. Don't check time since we might end up truncating
  //Check 0th dim separately in case of slicing
  if(my_data_in.get_dims(0) != space_range[1]-space_range[0]) err += 1;
  for(size_t i=1; i< my_data_in.get_dims()-1; i++)  if(my_data_in.get_dims(i) != source_sizes[i]) err += 1;
  if(err > 0){
    my_error_print("Source sizes do not match given array. Aborting read", mpi_info.rank);
    return 1;
  }

  bool simple_slice = false;
  int source_advance = 1;
  
  //Simple slices are those where we're slicing only the last spatial dimension, so we can just stop the read early. Since we only slice in x, that means either one and only one spatial dim, or no slicing.
  simple_slice = (my_data_in.get_dims() ==2 || source_sizes[0] == space_range[1]-space_range[0]);

  //Get SPACE advance for source, i.e advance between one time read and the next
  for(size_t i=0; i< grid_n_dims; i++){
    source_advance *= source_sizes[i];
  }

  //If we have to flatten, we need somewhere to put the flat data
  my_type * flat_data = nullptr;
  if(flatten_on >=0){
    flat_data = (my_type*) malloc(sizeof(my_type)*source_advance);
  }

  size_t i;
  size_t last_report=0;
  size_t report_interval = (time_range[1]-time_range[0])/10;
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;
    //Say we want to report 10 times over the list, or every 20th file if  more than 200.

  my_data_in.space[0] = space_range[0];
  my_data_in.space[1] = space_range[1];
  
  size_t rows=0; //Accumulated blocks only: number of rows to read in current file
  size_t total_reads=0; //Current number of reads done

  //Now loop over files getting data and time info
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
        //Complicated slicing, this is UNTESTED
        my_data_in.populate_complex_slice(my_ptr, 1, val, source_sizes);
      }
      //One time has been read
      total_reads++;
    }
    else{
      //Read all rows in file, but not more than time[2] rows total
      rows = block->dims[block->ndims-1];
      if(total_reads + rows >= time_range[2]) rows = time_range[2]- total_reads;

      read_acc_time(my_data_in, handle, total_reads, rows);
      size_t val[1] = {total_reads};
      //Handle each row in turn for simplicity
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
        //Update row to fill
        val[0] ++;
        //Update position in data file by correct amount
        my_ptr += source_advance;
      }
      //"Rows" times have been read
      total_reads+= rows;
    }
    sdf_close(handle);
    //Stop after time_range[2] total reads for accumulated data
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
    //Trim array to match number of lines actually read
    my_data_in.resize(1, total_reads);
    
    //Set last time, note we just resized the axis
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


bool reader::read_distrib(data_array & my_data_in, std::string dist_id, int dump_number){
/** \brief Read distribution function
*
*Reads the distribution function dist_id into data_in. ID has dist_fn removed and is trimmed to 10 chars max. If data type does not match compiled type we convert, assuming it is either a double or float type. my_data_in should be an array of correct size
*/
  if(!my_data_in.is_good()){
    my_error_print("Cannot read into invalid array", mpi_info.rank);
    return 1;
  }
  sdf_file_t * handle;
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

  //Copy the data, converting if necessary
  size_t total_els = my_data_in.get_total_elements();
  if(block->datatype != my_sdf_type) my_data_in.populate_data((other_type*)block->data, total_els);
  else my_data_in.populate_data((my_type*)block->data, total_els);
  
  //Copy the other info
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
    //Copy axes, converting if necessary
    if(block->datatype != my_sdf_type) std::copy((other_type *) block->grids[i], (other_type *) block->grids[i] + len, ax_ptr);
    else std::copy((my_type *) block->grids[i], (my_type *) block->grids[i] + len, ax_ptr);
  }
  
  return 0;

}

