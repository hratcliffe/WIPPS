//
//  data_array.cpp
//  
//
//  Created by Heather Ratcliffe on 03/08/2016.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <complex.h>
#include <limits.h>
#include <cmath>
#include <cstring>
#include "support.h"
#include "main_support.h"
#include "my_array.h"
#include "data_array.h"


extern const mpi_info_struc mpi_info;

void data_array::construct(){
/** \brief Common constructor logic
*
*Sets fields for empty, dimensionless array
*/
  axes=nullptr;
  time[0]=0; time[1]=1;
  space[0]=0; space[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  B_ref = 0;
  
}

void data_array::alloc_ax(const size_t els){
/* \brief Allocate axis memory
*
* Alocate memory for axes
*/
  if(els > 0 && els <= this->n_dims*MAX_SIZE){
    axes=(my_type*)calloc((els),sizeof(my_type));
  }else{
    if(els != 0) my_print("Array size exceeds max. Axes alloc failed", mpi_info.rank);
    else my_print("Array size cannot be 0. Axes alloc failed", mpi_info.rank);
  }
}

data_array::data_array() : my_array(){
/** \brief Default constructor
*
*Create an empty, dimensionless array
*/
  construct();
}

data_array::data_array(size_t nx, size_t ny, size_t nz, size_t nt) : my_array(nx, ny, nz, nt){
/**\brief Construct 1-4 d array
*
*Adds axes to a normal rectangular my array of correct size
*/
  construct();
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
}

data_array::data_array(size_t n_dims, size_t * dims) : my_array(n_dims, dims){
/** \brief Construct arbitrary dimension array
*
*Construct an array of n_dims rank with dimensions dims
*/
  construct();
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
}

data_array::data_array(std::string filename, bool no_version_check){
/**\brief Create data array from file
*
* Create a data array by reading from the named file. If the file does not exist no memory is allocated. Otherwise it reads the dimensions and sets itself up accordingly
*/

  std::fstream infile;
  infile.open(filename, std::ios::in|std::ios::binary);
  
  construct();
  if(!infile.is_open()){
    my_print("File open or access error");
    return;
  }
  std::vector<size_t> dims_vec = my_array::read_dims_from_file(infile, no_version_check);
  size_t n_dims_in = dims_vec.size();
  //Now we have the dimensions, construct
  size_t total_data=1, total_axes=0;
  if(n_dims_in >0){
    
    this->n_dims = n_dims_in;
    this->dims = (size_t*)malloc(n_dims*sizeof(size_t));
    for(size_t i=0;i<n_dims;i++){
      dims[i] = dims_vec[i];
      total_axes += dims[i];
      total_data *= dims[i];
      
    }
    
    data=(my_type*)calloc(total_data,sizeof(my_type));
    axes=(my_type*)calloc(total_axes,sizeof(my_type));

    //Finally read in data and axes using normal routines
    infile.seekg(0, std::ios::beg);
    bool err= this->read_from_file(infile, no_version_check);
    if(err){
      my_print("IO error, could not read", mpi_info.rank);
      if(axes) delete axes;
      if(data) delete data;
    }
  
  }else{
    my_print("Invalid dimensionality in input file", mpi_info.rank);
  }
}

data_array::~data_array(){
/**\brief Destructor 
*
*Free axis memory
**/
  if(axes) free(axes);
  axes = nullptr; // technically unnecessary as desctructor deletes members.
}

data_array & data_array::operator=(const data_array& src){
/** \brief Copy assignment
*
*Set this equal to src by copying src including data
*/
  if(this->axes) free(axes);

  this->construct();
  my_array::operator=(src);

  if(this->dims){
    size_t els= this->get_total_axis_elements();
    alloc_ax(els);
    //Allocate axis memory
    if(axes) std::copy(src.axes, src.axes+els, this->axes);
    //Copy axes
    copy_ids(src);
  }
  return *this;

}

data_array::data_array(data_array && src) : my_array(src){
/** \brief Move constructor
*
*Move src to new location. Copies data pointers but does not reallocate memory
*/
  this->axes = src.axes;
  src.axes = nullptr;
  copy_ids(src);
}

data_array::data_array(const data_array &src) : my_array(src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data
*/
  construct();
  //Basic construction of additionals, already called base class copy constructor
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
  //Allocate axis memory
  if(axes) std::copy(src.axes, src.axes+els, this->axes);
  //Copy axes
  copy_ids(src);
}

my_type * data_array::disown_axes(){
/** \brief Disown and return axes pointer
*
*Surrenders ownership of memory pointed to by axes NB if this pointer is not kept an manually freed, memory will leak NB calling disown_data on a data_array will not orphan the axes, but will make them inaccessible using getter/setter functions
*/
  my_type * data_pointer = this->axes;
  this->axes = nullptr;
  return data_pointer;
}


void data_array::clone_empty(const data_array &src){
/** \brief Initialise this to match sizes of src
*
* Allocates empty memory but does not copy data or axes
*/
  if(axes) free(axes);
  construct();
  my_array::clone_empty(src);
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);

}

my_type * data_array::get_axis(size_t dim, size_t & length){
/**  \brief Get pointer to axis
*
*Returns pointer to given axis and its length. If axes don't exist or dimension is out of range, returns nullptr
*/

  if(!axes || (dim >= n_dims)) return nullptr;

  long index = get_axis_index(dim, 0);
  //Get index of 0th element
  length = get_dims(dim);
  if(index != -1) return axes + index;
  else return nullptr;

}

long data_array::get_axis_index(size_t dim, size_t pt)const{
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. Let this function do all bounds checks.
*/

  if(dim >=n_dims || pt >=get_dims(dim)) return -1;
  //Out of range error
  
  long offset = 0;
  for(size_t i=0; i< dim; i++) offset +=dims[i];
  return offset + pt;

}

my_type data_array::get_axis_element(size_t dim, size_t pt)const{
/** \brief Get axis value
*
*Returns value at pt in dimension dim if in range, else 0.0
*/

  long ind = get_axis_index(dim, pt);
  if(ind  != -1){
    return axes[ind];
  }else{
    return 0.0;
  }

}

bool data_array::set_axis_element(size_t dim, size_t pt, my_type val){
/** \brief Sets array element
*
*Sets elements at pt on dimension dim, @return 1 if out of range, 0 else.
*/

  long index = get_axis_index(dim, pt);
  if(index >= 0){
    axes[index] = val;
    return 0;
  }else{
    return 1;
  }

}

bool data_array::populate_axis(size_t dim, my_type * dat_in, size_t n_tot){
/** \brief Fill axis
*
*Populates axis with the total number of elements specified, as long as n_tot is less than the specified dim. Parameter needed so we can't overflow dat_in. @return 0 (sucess) 1 (error)
*/

  if(dim >=n_dims) return 1;
  //Out of range error

  size_t tot_els = dims[dim];
  if(n_tot < tot_els) tot_els = n_tot;
  //Min of n_tot and tot_els

  void * tmp = (void *) this->axes;
  if(!tmp) return 1;

  memcpy (tmp , dat_in, tot_els*sizeof(my_type));

  return 0;

}   
float data_array::get_res(size_t i){
/**Return resolution of axis on dimension i. Assumes linear etc etc. If axis is undefined or zero or one in length, return 1.0 */
  size_t len;
  my_type * axis = this->get_axis(i, len);
  if(axis && len >1) return std::abs(axis[0]-axis[dims[i]-1])/(float)(dims[i]-1);
  else return 1.0;

}

size_t data_array::get_total_axis_elements(){
/** \brief Return total axes length
*
*Sums number of total elements in all axes
*/
  size_t tot_els=0;

  for(size_t i=0; i<n_dims;++i) tot_els +=dims[i];

  return tot_els;

}

void data_array::make_linear_axis(size_t dim, float res, long offset){
/**\brief Make an axis
*
*Generates a linear axis for dimension dim, with resolution res, starting at value of  - offset*res @param dim Dimension to build axis for @param res Axis resolution @param offset Number of grid cells to shift downwards (leftwards) by
*/

  size_t len;
  my_type * ax_ptr = get_axis(dim, len);

  for(size_t i=0; i<len; i++){
    *(ax_ptr +i) = ((float) ((long)i-offset)) * res;
  }
}

bool data_array::write_to_file(std::fstream &file, bool close_file){
/** \brief Write data array to file
*
*First write the my_array section, see my_array::write_to_file then add the following
* Next_block axes
* Next_block Block_id Prev_block

IMPORTANT: the VERSION specifier links output files to code. If modifying output or order commit and clean build before using.
*/

  if(!file.is_open()) return 1;
  bool write_err=0;

  my_array::write_to_file(file);
  //call base class method to write that data.

  size_t tot_ax = get_total_axis_elements();

  size_t next_location = (size_t) file.tellg() + sizeof(size_t) + tot_ax*sizeof(my_type);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section

  file.write((char *) axes , sizeof(my_type)*tot_ax);
  //Add axes.
  if((size_t)file.tellg() != next_location) write_err=1;

  size_t ftr_start = next_location;

  next_location += sizeof(my_type)*3 + sizeof(size_t)*2 + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write((char *)time, sizeof(my_type)*2);
  file.write((char *)space, sizeof(size_t)*2);
  file.write((char *) &B_ref, sizeof(my_type));

  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  if(close_file) file.write((char*) & ftr_start, sizeof(size_t));
  //Finish with position of start of footer!
  return 0;

}

bool data_array::write_section_to_file(std::fstream &file, std::vector<my_type> limits, bool close_file){
/** \brief Write section to file
*
*Write section between given AXIS values to file. To use one dimension entire supply values less/greater than min and max axis values.
*/
  //Identify limits of segment from axes
  std::vector<size_t> index_limits = this->get_bounds(limits);
  return write_raw_section_to_file(file, index_limits, close_file);

}

bool data_array::write_raw_section_to_file(std::fstream &file, std::vector<size_t> index_limits, bool close_file){
/** \brief Print section of array to file
*
*Prints the section defined by the vector limits to supplied file. If limits are out of range the start and end are used. See data_array::write_to_file for format
*/

  if(!file.is_open()) return 1;
  bool write_err = 0;

  if(index_limits.size() != 2*n_dims){
    my_print("Limits vector size does not match array!", mpi_info.rank);
    return 1;
  }

  for(size_t i=0; i< index_limits.size(); i+=2){
    if(index_limits[i] >= dims[i/2]) index_limits[i] = 0;
    if(index_limits[i+1] >= dims[i/2]) index_limits[i+1] = dims[i/2];
    //Make sure are within ranges.
  }

  my_array::write_section_to_file(file, index_limits);
  //call base class method to write that data.

  size_t tot_ax = 0;
  for(size_t i=0; i< n_dims; i++){
    tot_ax +=(index_limits[2*i +1]-index_limits[2*i]);
  }
  
  size_t next_location = (size_t) file.tellg() + sizeof(size_t) + tot_ax*sizeof(my_type);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section

  size_t len;
  for(size_t i=0; i< n_dims; i++){
    file.write((char *) (get_axis(i, len)+index_limits[2*i]), sizeof(my_type)*(index_limits[2*i +1]-index_limits[2*i]));

  }
  //Add axes.

  if((size_t)file.tellg() != next_location) write_err=1;

  size_t ftr_start = next_location;
  //Start of ftr means where to start reading block, i.e. location of the next_location tag

  next_location += sizeof(my_type)*3 + sizeof(size_t)*2 + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write((char *) time, sizeof(my_type)*2);
  file.write((char *) space, sizeof(size_t)*2);
  file.write((char *) &B_ref, sizeof(my_type));

  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  if(close_file) file.write((char*) & ftr_start, sizeof(size_t));

  return 0;
}

bool data_array::write_closer(std::fstream &file){

  if(!file.is_open()) return 1;
  bool write_err = 0;

  size_t ftr_start = (size_t) file.tellg();
  //Start of ftr means where to start reading block, i.e. location of the next_location tag
  size_t next_location = ftr_start+ sizeof(char)*ID_SIZE +sizeof(size_t);

  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section, i.e. of file end val
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & ftr_start, sizeof(size_t));

  return write_err;
}

std::vector<size_t> data_array::get_bounds(std::vector<my_type> limits){

  std::vector<size_t> index_limits;

  if(limits.size() != 2*n_dims){
    my_print("Limits vector size does not match array!", mpi_info.rank);
    return index_limits;
  }

  index_limits.resize(n_dims*2);

  size_t len;
//  long index;
  my_type * ax_start;
  long where_val;
  for(size_t i=0; i< n_dims; i++){
    ax_start = get_axis(i, len);
    where_val = where(ax_start, len, limits[2*i]);
    if(where_val != -1) index_limits[i*2] = where_val;
    else index_limits[i*2] = 0;
    where_val = where(ax_start, len, limits[2*i + 1]);
    if(where_val != -1) index_limits[i*2 + 1] = where_val;
    else index_limits[i*2 + 1] = this->dims[i];
  }

  return index_limits;
}

bool data_array::read_from_file(std::fstream &file, bool no_version_check){
/** \brief Read data array file dump
*
*Read a file into a pre-sized data array See also data_array::data_array(std::string filename, bool no_version_check)
*/
/*
*First write the my_array section, see my_array::write_to_file then add the following
* Next_block axes
* Next_block Block_id Prev_block
*/

  bool err=1;

  if(file.good()) err=my_array::read_from_file(file, no_version_check);
  if(err){
    my_print("File read failed", mpi_info.rank);
    return err;
  }
  //call parent class to read data, checking we read id ok first

  size_t next_block=0; //, end_block=0, end_pos=0;
  file.read((char*) &next_block, sizeof(size_t));

  size_t tot_els = get_total_axis_elements();

  if(!err) file.read((char *) this->axes , sizeof(my_type)*(tot_els));
  //If we managed to get dims etc, continue on to get axes

  file.read((char*) &next_block, sizeof(size_t));

  //Position of next section
  file.read((char *) time, sizeof(my_type)*2);
  file.read((char *) space, sizeof(size_t)*2);
  file.read((char *) &B_ref, sizeof(my_type));

  file.read((char*) &next_block, sizeof(size_t));

  //Now read the block ID
  char id_in[ID_SIZE];
//  file.seekg(end_block+sizeof(size_t));
  if(file){
    file.read(id_in, sizeof(char)*ID_SIZE);
    strncpy(this->block_id, id_in, ID_SIZE);
  }
  //file.seekg(end_pos);
  return err;

}

bool data_array::fft_me(data_array & data_out){
#ifndef NO_FFT
/** \brief FFT data_array
*
* Data and axes in this object are FFT'd using FFTW and stored into the instance pointed to by data_out. Data_out must be created with correct dimensions first, but we check and return error (1) if it is not so. \todo extract and make friend?
*/

  if(!data_out.is_good()){
    my_print("Output array for FFT undefined", mpi_info.rank);
    return 1;
  }
  if(data_out.n_dims != this->n_dims){
    my_print("Wrong output dimensions for FFT", mpi_info.rank);
    return 1;
  }
  for(size_t i=0; i<n_dims;++i){
    if(data_out.dims[i] != this->dims[i]){
      my_print("Wrong output dimensions for FFT", mpi_info.rank);
      return 1;
    }
  }

  data_out.copy_ids(*this);
  size_t total_size=1; /* Total number of elements in array*/
  for(size_t i=0; i<n_dims;++i) total_size *= dims[i];

  ADD_FFTW(plan) p;
  cplx_type *out;
  my_type * in, *result;

  size_t output_size=this->get_total_elements()/dims[0]*(dims[0]/2+1);
  //Size of r2c transform output
  //This has been checked manually with valgrind and is CORRECT

  in = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);
  //my_type should match the used FFTW library, so no type conversion necessary
  out = (cplx_type *) ADD_FFTW(malloc)(sizeof(cplx_type) * output_size);

  result = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * output_size);

  int *rev_dims;
  rev_dims = (int *) malloc(n_dims*sizeof(size_t));
  for(size_t i=0; i< n_dims; i++) rev_dims[i] =dims[n_dims-1-i];
  //Construct dims array in reverse order as we're using the opposite majority to FFTW
  p = ADD_FFTW(plan_dft_r2c)((int)n_dims, rev_dims, in, out, FFTW_ESTIMATE);
  free(rev_dims);
/*
  if(n_dims == 1 || (n_dims == 2 && dims[1] == 1) ){
    p = ADD_FFTW(plan_dft_r2c_1d)(dims[0], in, out, FFTW_ESTIMATE);

  }else if(n_dims == 2){
    p = ADD_FFTW(plan_dft_r2c_2d)(dims[1], dims[0], in, out, FFTW_ESTIMATE);

  }else if(n_dims == 3){
    p = ADD_FFTW(plan_dft_r2c_3d)(dims[2], dims[1], dims[0], in, out, FFTW_ESTIMATE);

  }else{
    my_print("FFT of more than 2-d arrays not added yet", mpi_info.rank);

    return 1;
  }*/

  //copy data into in. Because the plan creation changes in, so we don't want to feed our actual data array in, and it's safer to run the plan with the memory block it was created with
  std::copy(this->data, this->data+total_size, in);

  ADD_FFTW(execute)(p);
  //Execute the plan

  cplx_type * addr;
  addr = out;
  //because double indirection is messy and cplx type is currently a 2-element array of floats

  for(size_t i=0; i< output_size; i++){
    *(result+i) = (my_type)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1]));
    addr++;
  }
  //Absolute square of out array to produce final result of type my_type

  ADD_FFTW(free)(in);
  ADD_FFTW(free)(out);
  //Free as early as possible...
  
  bool err=false;
  err = data_out.populate_mirror_fastest(result, total_size);
  //Copy result into out array

  if(err){
    my_print("Error populating result array", mpi_info.rank);
    return 1;
  }

  ADD_FFTW(destroy_plan)(p);
  ADD_FFTW(free)(result);
  //Destroy stuff we don't need

  my_type * tmp_axis;
  float N2, res;
  size_t len;

  for(size_t i=0;i<n_dims;++i){
  //Loop over the dimensions and construct each axis in place. We KNOW now that data_out has correct dimensions. We checked before getting here. We don't need to check len. It is the same as dims[i] each time. Perhaps we will anyway? :)
    tmp_axis = data_out.get_axis(i, len);
    if(len != dims[i]) return 1;

    N2 = ((float) dims[i])/2.0;
    res = this->get_res(i);
    for(size_t j= 0; j< dims[i]; j++) *(tmp_axis + j) = pi * ((float)j -N2)/N2/res;
  }

  return 0;

#endif
return 1;
//Return for noFFT version
}

bool data_array::populate_mirror_fastest(my_type * result_in, size_t total_els){
/** \brief Copy FFTW data into array

*
*For real data an FFT has a redundant half so FFTW returns array od size (dims[0]/2 +1)*dims[...]*dims[n-1] (Note that our first dim is FFTWs last). Since we want all k we have to mirror this to obey H(-f, -g) = H(f, g). For simplicity we enforce this on k_x, omega pair only. Data is assumed unshifted. Should work for 1-3 dimensions. 1st dimension will be returned shifted to 0-in-centre \todo Which side is negative k?
*/
  size_t last_size = floor(dims[0]/2) + 1;
  size_t num_strides = total_els/dims[0];//Always integer

  size_t num_seconds = this->n_dims >= 2? dims[1]: 1;//Size of "second dimension"
  size_t slice_sz = num_seconds*dims[0];

memset(data, 0, total_els);

  for(size_t i=0; i< num_strides; i++){
    //Reverse the result into place (left side of a 2-d)
    std::reverse_copy(result_in+ i*last_size, result_in+(i+1)*last_size-1, data+i*dims[0]+1);
  }

  //Now fill the redundant 2-d side. Zero stays as zero, we flip the rest. This is to give us "correct" 2-D result when we take only +ve omega
  for(size_t j=0; j< num_strides/num_seconds; j++){
    std::reverse_copy(data+1+j*slice_sz, data+last_size+j*slice_sz, data+last_size-1+j*slice_sz);
    for(size_t i=1; i< num_seconds; i++){
      std::reverse_copy(data+i*dims[0]+1+j*slice_sz, data+i*dims[0]+last_size+j*slice_sz, data+j*slice_sz+(num_seconds-i)*dims[0]+last_size-1);
    }
  }
  
  //do shifts for all but 0th dim, this is already done by steps above
  size_t shft = 0;
  for(size_t i=1; i< n_dims; i++){
    shft = this->get_dims(i)/2;
    this->shift(i, shft, 0);
  }

  
  if(n_dims >=3){
    //Reverse omega vs k_x on -ve k_x side only.
    for(size_t i=0; i<dims[0]; i++){
      for(size_t j=0; j<dims[1]; j++){
        for(size_t k=0; k<dims[2]/2; k++){
        
          *(data+ i + j*dims[0]+(dims[2]-1-k)*dims[1]*dims[0]) = *(data+ i + j*dims[0]+ k*dims[1]*dims[0]);
        }
      }
    }
  
  }
  
  return 0;
}

void data_array::copy_ids( const data_array &src){
/** Copies ID fields from src array to this */

  strncpy(this->block_id, src.block_id, ID_SIZE);
  //Strncpy for extra safety
  for(size_t i=0; i < 2; ++i) this->time[i] = src.time[i];
  //  std::copy(src->time, src->time + 2, this->time);
  for(size_t i=0; i < 2; ++i) this->space[i] = src.space[i];
  this->B_ref = src.B_ref;
}

bool data_array::check_ids( const data_array & src){
/** Checks ID fields match src */

  bool err=false;
  if(strcmp(this->block_id, src.block_id) != 0) err =true;
  for(size_t i=0; i< 2; i++) if(src.time[i] != this->time[i]) err=true;
  for(size_t i=0; i < 2; ++i) if(this->space[i] != src.space[i]) err=true;
  err |= this->B_ref != src.B_ref;
  return err;
}

bool data_array::resize(size_t dim, size_t sz, bool verbose){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Similarly for axis elements. See my_array::resize() for more.
*/

  size_t old_sz = this->get_dims(dim);
  //call my_array::resize to resize data...
  bool err = my_array::resize(dim, sz, verbose);

  if(!err){
    //success! Do axes!
    my_type * new_ax;
    size_t part_sz = 0;

    if(dim == n_dims-1){
      //special case as we can shrink and maybe grow without copy
      for(size_t i=0; i<dim; ++i) part_sz+= dims[i];
      //product of all other dims
      new_ax = (my_type *) realloc((void*) this->axes, (part_sz+sz)*sizeof(my_type));
      if(!new_ax){
        my_print("Failed to reallocate axes", mpi_info.rank);
        return 1;
        //failure. leave as was.
      }
      long new_els = (sz - dims[dim]);

      if(new_els > 0) memset((void*)(new_ax + part_sz), 0, new_els*sizeof(my_type));
      //zero new elements
      axes = new_ax;
    }else{
      //have to allocate a new block and copy across.
      for(size_t i=0; i<dim; ++i) part_sz+= dims[i];
      for(size_t i=dim+1; i<n_dims; ++i) part_sz+= dims[i];

      new_ax=(my_type*)calloc(part_sz+sz,sizeof(my_type));
      size_t els_to_copy = 0, old_starts=0, new_starts=0;
      for(size_t i=0;i<n_dims; ++i){
        els_to_copy = dims[i];
        if(i == dim && sz< old_sz) els_to_copy = sz;
        std::copy(axes+old_starts, axes+old_starts+els_to_copy, new_ax+new_starts);

        if(i!= dim) old_starts += dims[i];
        else old_starts +=old_sz;
        new_starts +=els_to_copy;
      
      }
      free(axes);
      axes = new_ax;

    }
    return 0;
  }else{
  
    return err;
  
  }


}

bool data_array::shift(size_t dim, long n_els, bool axis){
/** Shift array on dim dim by n_els
*
*@param axis Whether to shift the corresponding axis
*/
  if(dim >= n_dims) return 0;
  bool err;
  err = my_array::shift(dim, n_els);

  if(axis){
    long sign_n = (n_els<0? -1: 1);
    n_els = sign_n >0? (std::abs(n_els)%dims[dim]):dims[dim]-(std::abs(n_els)%dims[dim]) ;
  //By now n_els guaranteed >= 0

    size_t len=0;
    my_type * ax = get_axis(dim, len);
    my_type * new_data=(my_type*)malloc(len*sizeof(my_type));
    if(!new_data){
      my_print("Error allocating spare memory for shift", mpi_info.rank);
      return 1;
    
    }
    long actual_shift = 0;
    //Rotate the axis
    std::copy(ax,ax+len, new_data);
    for(size_t j=0; j< len; j++){
      //This is inefficient but is minimmaly changed from shift for whole array
      actual_shift = (j+n_els >= dims[dim])? (j+n_els-dims[dim]): j+n_els;
      actual_shift += (actual_shift < 0 ? dims[dim]:0);
      std::copy(new_data+j, new_data+(j+1), ax+actual_shift);
    }
    if(new_data) free(new_data);
  }
  return err;

}

data_array data_array::total(size_t dim){
/** \brief Total array on dim dim
*
*Returns a new array of rank this->n_dims -1, containing data summed over entire range of dimension dim. See also data_array data_array::total(size_t dim, my_type min, my_type max)
*/
  //Wraps total over range, passing limits outside axis ranges
  my_type min, max;
  min = this->get_axis_element(dim, 0) - 1;
  max = this->get_axis_element(dim, this->get_dims(dim)-1) + 1;
  return this->total(dim, min, max);

}

data_array data_array::total(size_t dim, my_type min, my_type max){

  size_t len, min_ind = 0, max_ind = this->get_dims(dim)-1;
  int where_val;
  my_type * ax_start = get_axis(dim, len);
  where_val = where(ax_start, len, min);
  if(where_val != -1) min_ind = where_val;
  where_val = where(ax_start, len, max);
  if(where_val != -1) max_ind = where_val;


  return total(dim, min_ind, max_ind);

}

data_array data_array::total(size_t dim, size_t min_ind, size_t max_ind){

/** \brief Total along dim dim
*
*Returns a new array of rank this->n_dims -1, containing data summed over between axis values of min and max on dimension dim. See also data_array data_array::total(size_t dim)
*/
  
  if(dim >= this->n_dims){
    return data_array();
  }
  
  size_t * new_dims;
  new_dims = (size_t *) calloc((this->n_dims-1), sizeof(size_t));
  for(size_t i=0, i2=0; i< this->n_dims; i++, i2++){
    if(i == dim){
      i2--;
      continue;
    }
    new_dims[i2] = dims[i];
  }

  data_array new_array = data_array(this->n_dims-1, new_dims);
  
  flatten_fortran_slice(this->data, new_array.data, this->n_dims, this->dims, dim, min_ind, max_ind);
  
  //Copy axes
  my_type * ax_new;
  size_t len2, len;
  for(size_t i=0, i2=0; i2< n_dims -1; i++, i2++ ){
    if(i == dim) i++;
    ax_new = new_array.get_axis(i2, len);
    std::copy(this->get_axis(i, len2),this->get_axis(i, len2)+len, ax_new);
  }
  new_array.copy_ids(*this);
  free(new_dims);
  return new_array;
}

data_array data_array::average(size_t dim){

  size_t old_dim = this->get_dims(dim);
  data_array new_arr = this->total(dim);
  new_arr.divide(old_dim);
  return new_arr;
}


data_array data_array::average(size_t dim, my_type min, my_type max){
  size_t len, min_ind=0, max_ind = this->get_dims(dim)-1;
  int where_val;
  my_type * ax_start = get_axis(dim, len);
  where_val = where(ax_start, len, min);
  if(where_val != -1) min_ind = where_val;
  where_val = where(ax_start, len, max);
  if(where_val != -1) max_ind = where_val;

  return average(dim, min_ind, max_ind);
}

data_array data_array::average(size_t dim, size_t min_ind, size_t max_ind){
  data_array new_arr = total(dim, min_ind, max_ind);
  new_arr.divide(max_ind - min_ind);
  return new_arr;
}

bool data_array::subtract(const data_array& rhs){
  
  if(this->get_dims() != rhs.get_dims()) return 1;
  for(size_t i=0; i< n_dims; i++) if(this->get_dims(i) != rhs.get_dims(i)) return 1;

/** \todo Why is this copying the axes????*/
  size_t els= get_total_axis_elements();
  if(this->axes) std::copy(rhs.axes, rhs.axes+els, this->axes);
    //Copy axes

  els = get_total_elements();
  for(size_t i=0; i< els; i++) *(this->data+i) -= rhs.get_element_from_index(i);
  return 0;
}

