//
//  my_array.cpp
//  
//
//  Created by Heather Ratcliffe on 21/09/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <cstring>
#include <numeric>
#include <algorithm>
#include "support.h"
#include "my_array.h"
#include "data_array.h"

/********Basic setup and allocation functions ****/

my_array::my_array(){
/** Default constructor*/
  construct();
}

void my_array::construct(){
/** \brief Shared contructor code
*
*Sets default values
*/
  n_dims = 0;
  data = nullptr;
  dims = nullptr;

}

void my_array::alloc_all(const size_t n_dims, const size_t * const dims){
/** \brief Takes care of memory allocation
*
* Allocate dims array and data. If any size is zero, any exceeds MAX_SIZE, or overall exceeds MAX_SIZE_TOT, print error and stop. NB inputs will not be modified, will be copied
@param n_dims Rank of array to allocate
@param dims Dimensions of array to allocate
*/

  size_t tot_dims = 1;
  bool too_large = false;
  for(size_t i=0; i< n_dims; i++){
    tot_dims *= dims[i];
    if(dims[i] > MAX_SIZE) too_large = true;
  }
  //tot_dims now 0 if any dim is 0
  if(tot_dims==0){
    my_error_print("Array cannot have 0 dim", mpi_info.rank);
    return;
  }
  if(too_large || tot_dims > MAX_SIZE_TOT){
    if(too_large) my_error_print("Array size exceeds max per-dimension size of "+mk_str(MAX_SIZE), mpi_info.rank);
    if(tot_dims> MAX_SIZE_TOT) my_error_print("Array size exceeds max overall size of "+mk_str(MAX_SIZE_TOT), mpi_info.rank);
    return;
  }
  
  this->n_dims = n_dims;
  this->dims = (size_t*)malloc(n_dims*sizeof(size_t));
  if(dims){
  //If allocation succeeds, we can copy in dims data and allocate space for data, else stop
    std::copy(dims, dims+n_dims, this->dims);
    data = (my_type*) calloc(tot_dims, sizeof(my_type));
  }

}

my_array::my_array(size_t nx, size_t ny, size_t nz, size_t nt){
/** \brief 1 to 4 d rectangular array creator
*
*Sets up a n-d rectangular array for n = 1 to 4. Helper avoids user having to construct size_t array of dims. Default values mean any dimensions not supplied will be 0
@param nx Size of x dimension
@param ny Size of y dimension, default 0
@param nz Size of z dimension, default 0
@param nt Size of t dimension, default 0
*/
  construct();
  
  size_t * dims_in;
  size_t n_dims_act = 4;

  //Construct the dims array to pass to normal constructor
  dims_in = (size_t*) malloc(n_dims_act*sizeof(size_t));
  dims_in[0] = nx;
  dims_in[1] = ny;
  dims_in[2] = nz;
  dims_in[3] = nt;
  
  //Find highest non-zero dimension
  for(size_t i=4; i>0; i--){
    if(dims_in[i-1] == 0) n_dims_act --;
    else break;
  }
  //Allocate as normal using dims_in[0:n_dims_act]
  alloc_all(n_dims_act, dims_in);
  free(dims_in);
}

my_array::my_array(size_t n_dims, size_t * dims ){
/** \brief Arbitrary dim rectangular array
*
*Sets up internals of array including memory allocation
@param n_dims Rank of array to create
@param dims Array of dimensions of array to create
*/
  construct();
  alloc_all(n_dims, dims);
}

my_array::~my_array(){
/** \brief Destructor
*
*Clean up explicit allocations
*/
  if(data) free(data);
  data=nullptr;
  if(dims) free(dims);
  dims=nullptr;

}

/********Technical stuff making my_array a proper "object" ****/

my_array & my_array::operator=(const my_array& src){
/** \brief Copy assignment
*
*Sets this equal to a (deep) copy of source
@param src Array to copy from
@return Copy of input
*/

  //Trap self-assigning or bad copy before destructing
  if(&src == this || !src.is_good()) return *this;
  
  if(this->data) free(data);
  if(this->dims) free(dims);
  //Clean up in case this was already allocated
  my_array::construct();
  
  //Construct this and copy all necessary
  alloc_all(src.n_dims, src.dims);
  if(this->dims){
    size_t tot_els = this->get_total_elements();
    if(this->data && src.data) std::copy(src.data, src.data + tot_els, this->data);
  }
  return *this;

}

my_array::my_array(my_array && src){
/** \brief Move constructor
*
*Move src to a new instance i.e. copy fields but don't move memory. Src becomes empty afterwards
@param src Array to move
*/

  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims
  
  this->n_dims = src.n_dims;
  src.n_dims = 0;
  //Inherit memory and unlink from src
  this->dims = src.dims;
  src.dims = nullptr;
  this->data = src.data;
  src.data = nullptr;
}

my_array::my_array(const my_array &src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data
@param src Array to copy from
*/

  construct();
  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims
  
  alloc_all(src.n_dims, src.dims);
  size_t tot_els = this->get_total_elements();
  if(this->data && src.is_good()) std::copy(src.data, src.data + tot_els, this->data);

}

bool my_array::operator==(const my_array &rhs)const{
/** \brief Equality operator
*
* Check this is equal to rhs. Since copies are always deep, we check values, not data pointers
@param rhs Array to compare to
@return Boolean true if equal, false else
*/
  if(this->get_dims() != rhs.get_dims()) return false;
  for(size_t i=0; i< this->get_dims(); i++) if(this->get_dims(i) != rhs.get_dims(i)) return false;
  //Check each element in 1-D array
  for(size_t i=0; i< this->get_total_elements(); i++) if(*(this->data + i) != *(rhs.data + i)) return false;

  return true;
}

/********Helpers for working with my_array ****/

my_type * my_array::disown_data(){
/** \brief Disown and return data pointer
*
*Surrenders ownership of memory pointed to by data, nullifies dimensions and returns pointer. NB if this pointer is not kept and manually freed, memory will leak. This array will then be an empty array
@return Pointer to the disowned data array
*/
  my_type * data_pointer = this->data;
  this->data = nullptr;
  this->n_dims = 0;
  free(this->dims);
  this->dims = nullptr;
  return data_pointer;
}

void my_array::clone_empty(const my_array &src){
/** \brief Initialise this to same sizes as src
*
* This will be a valid empty array of size matching src.
@param src Array to copy dims from
*/

  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims

  //Delete any previous memory allocations and reset sizes
  if(data) free(data);
  if(dims) free(dims);
  construct();

  alloc_all(src.n_dims, src.dims);
  
}

bool my_array::copy_data(my_type * destination)const{
/** \brief Copy the data into destination array
*
*Data is not lost, but a direct copy is made. Destination size will NOT be checked, the entirety of data is copied.
@param[out] destination Pointer to destination to copy to
@return 0, error checking not yet implemented
*/
  std::copy(this->data, this->data+this->get_total_elements(), destination);
  return 0; //Add error checking??
}

void my_array::zero_data(){
/** \brief Reset data to 0
*
*
*/
  memset(data, 0.0, this->get_total_elements()*sizeof(my_type));
}

/********Indexers, getters and setters ****/

size_t my_array::get_dims()const{
/** \brief Return rank of array
*/
  return n_dims;
}
size_t my_array::get_dims(size_t dim)const{
/** \brief Return size of dimension dim
*
*/
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

long my_array::get_index(size_t n_dims, size_t * inds_in)const{
/** \brief Convert n-d index into 1-d index
*
*This requires passing integer array and loop so is slower than the dedicated functions (get_index(nx, ...)) that follow
*/
  //Check rank and bounds
  if(n_dims != this-> n_dims) return -1;
  for(size_t i=0; i< n_dims; i++){
      if(inds_in[i] >=dims[i]) return -1;
  }
  //Construct 1-d index from array
  long ind = inds_in[0];
  size_t subdims = dims[0];
  for(size_t i=1; i< n_dims; i++){
    ind += subdims*inds_in[i];
    subdims *= dims[i];
  }
  return ind;
  
}

long my_array::get_index(size_t nx)const{
/** \brief Get index of element at nx
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, the index into backing data array. This function is called often so we make it as simple as possible and write one for each number of args.
*/
/* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>
Etc */

#ifdef DEBUG_DIMS
  if(n_dims != 1){
  //Debug message
    my_error_print("Wrong array dimension, attempting 1 with "+mk_str(n_dims), mpi_info.rank);
    return -1;
  }
#endif
  if((nx < dims[0])){
    //Long should be big enough
    return (long) nx;
  }else{
#ifdef DEBUG_RANGE
    my_error_print("Array index out of range", mpi_info.rank);
    throw std::out_of_range ("1D array indexer");
#endif
    return -1;
  }

}
long my_array::get_index(size_t nx, size_t ny)const{
/**
  See my_array::get_index(size_t nx) const
*/
#ifdef DEBUG_DIMS
  if(n_dims != 2){
    my_error_print("Wrong array dimension, attempting 2 with "+mk_str(n_dims), mpi_info.rank);
    return -1;
  }
#endif
  if((nx < dims[0]) && (ny<dims[1])){
    return (long) (ny*dims[0] + nx);
  }else{
#ifdef DEBUG_RANGE
    my_error_print("Array index out of range", mpi_info.rank);
    throw std::out_of_range ("2D array indexer");
#endif
    return -1;
  }
}
long my_array::get_index(size_t nx, size_t ny, size_t nz)const{
/**
  See my_array::get_index(size_t nx) const
*/
#ifdef DEBUG_DIMS
  if(n_dims != 3){
    my_error_print("Wrong array dimension, attempting 3 with "+mk_str(n_dims), mpi_info.rank);
    return -1;
  }
#endif
  if((nx < dims[0]) && (ny<dims[1]) && ((nz<dims[2]))){
    return (nz*dims[1]+ ny)*dims[0] + nx;
  }else{
#ifdef DEBUG_RANGE
    my_error_print("Array index out of range", mpi_info.rank);
    throw std::out_of_range ("3D array indexer");
#endif
    return -1;
  }
}
long my_array::get_index(size_t nx, size_t ny, size_t nz, size_t nt)const{
/**
  See my_array::get_index(size_t nx) const
*/
#ifdef DEBUG_DIMS
  if(n_dims != 4){
    my_error_print("Wrong array dimension, attempting 4 with "+mk_str(n_dims), mpi_info.rank);
    return -1;
  }
#endif
  if((nx < dims[0]) && (ny<dims[1])&& ((nz<dims[2]))&& ((nt<dims[3]))){
    return ((nt*dims[2] +nz)*dims[1]+ ny)*dims[0] + nx;
  }else{
#ifdef DEBUG_RANGE
    my_error_print("Array index out of range", mpi_info.rank);
    throw std::out_of_range ("Arb-D array indexer");
#endif
    return -1;
  }
}

std::vector<size_t> my_array::get_indices_from_offset(size_t offset)const{
/** \brief Get vector of indices from 1-d index
*
*Takes a 1-d index and returns the n-d index vector. If offset is out of range, empty vector is returned.
*/
  std::vector<size_t> pos;
  if(offset > this->get_total_elements()) return pos;

  pos.resize(n_dims);
  size_t tmp_offset = offset;
  
  size_t stride = 1;
  for(size_t i=0; i< n_dims-1; i++) stride*=dims[i];

  for(size_t i=n_dims-1; i>0; i--){
    pos[i] = tmp_offset/stride;
    //INTEGER division!!!
    tmp_offset = tmp_offset%stride;
    stride/=dims[i-1];
  }
  pos[0] = tmp_offset/stride;
  return pos;
}

my_type my_array::get_element(size_t nx)const{
/** \brief Get element
*
* Return element at nx. Out of range etc will return 0.0*/
  long ind = get_index(nx);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t nx, size_t ny)const{
/** As my_array::get_element(size_t nx) const but for 2-D arrays*/
  long ind = get_index(nx, ny);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz)const{
/** As my_array::get_element(size_t nx) const but for 3-D arrays*/
  long ind = get_index(nx, ny, nz);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz, size_t nt)const{
/** As my_array::get_element(size_t nx) const but for 4-D arrays*/
  long ind = get_index(nx, ny, nz, nt);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t n_dims, size_t* dims)const{
/** As my_array::get_element(size_t nx) const but for arbitrary dimension arrays. Supply n_dims and an array of the required indices*/
  long ind = get_index(n_dims, dims);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element_from_index(size_t ind)const{
/** \brief Get element by 1-d offset
*
*Returns the element at index in the 1-D backing array. Ind should be found using one of the get_index options
*/
  if(ind < get_total_elements()){
    return *(data+ind);
  }else{
    return 0.0;
  }
}

size_t my_array::get_total_elements()const{
/** Return total size of array 
@return The total number of data elements in array*/
  size_t tot_els=1;

  for(size_t i=0; i<n_dims;++i) tot_els *=dims[i];
  return tot_els;

}

bool my_array::set_element(size_t nx, my_type val){
/** \brief Sets array element
*
*Sets elements at nx. @return 1 if nx is out of range or array is not rank 1, 0 else.
*/

  long index = get_index(nx);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
}
bool my_array::set_element(size_t nx, size_t ny, my_type val){
/**As my_array::set_element(size_t nx, my_type val) for 2-D arrays*/

  long index = get_index(nx, ny);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
}
bool my_array::set_element(size_t nx, size_t ny, size_t nz, my_type val){
/**As my_array::set_element(size_t nx, my_type val) for 3-D arrays*/

  long index = get_index(nx, ny, nz);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
}
bool my_array::set_element(size_t nx, size_t ny, size_t nz, size_t nt, my_type val){
/**As my_array::set_element(size_t nx, my_type val) for 4-D arrays*/

  long index = get_index(nx, ny, nz, nt);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
}
bool my_array::set_element(size_t n_dims_in, size_t * indices, my_type val){
/**As my_array::set_element(size_t nx, my_type val) for N-D arrays, using array of indexes*/

  long index = get_index(n_dims_in, indices);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
}

/********Data fillers, file IO ****/

template<typename T> bool my_array::populate_data(T dat_in, size_t n_tot){
/** \brief Fill array
*
*Populates this array from dat_in. n_tot should be the total number of elements in dat_in. The smaller of n_tot and the total number of elements in this array are copied. dat_in must match this array in row-column ordering and rank. 
@param dat_in Source of copy
@param n_tot Total number of elements to copy
@return 0 (sucess) 1 (error).
*/
  size_t tot_els = get_total_elements();
#ifdef DEBUG_DIMS
  if(tot_els != n_tot){
    my_error_print("Attempting fill with mismatched array"+mk_str(tot_els)+" vs "+mk_str(n_tot), mpi_info.rank);
    return 1;
  }
#endif
  if(n_tot < tot_els) tot_els = n_tot;
  //Use min of n_tot, tot_els
  std::copy(dat_in, dat_in+tot_els, this->data);
  return 0;

}
//Create only my_type and other_type versions of the populate function though
template bool my_array::populate_data(my_type *, size_t);
template bool my_array::populate_data(other_type *, size_t);

bool my_array::populate_slice(my_type * dat_in, size_t offsets_n_dims, size_t * offsets){
/** \brief Populate a slice of array from the input array.
*
* Fill an array slice, that is a section of rank m < n_dims, with some finite offset in dimensions from m to n_dims only. offsets_n_dims = n_dims - m is the size of the array of offsets provided. E.g. to fill a row of a 3-d array call with n_dims_in = 3-1=2 and offsets={column, plane}. Or to fill an array of shape (x, y, t) at a single time value t_0, use n_dims_in=1, offsets={t_0}.
 @param dat_in pointer to data 
 @param n_dims Dimensionality of input (must be less than dimension of array) 
 @param offsets Offsets in the other dimensions
 @return 0 for success, 1 for error
*/

  if(offsets_n_dims >= n_dims) return 1;

  long indx;
  //Infer size of dat_in as product of sizes of all dims not sliced on
  size_t sz_in = 1;
  size_t * inds_arr;
  inds_arr = (size_t *) calloc(n_dims,sizeof(size_t));
  for(size_t i=0; i<n_dims-offsets_n_dims; i++){
    sz_in *=dims[i];
  }
  //Location to copy to is 0 in all the matched dimensions and matches offsets in all the rest
  for(size_t i = n_dims-offsets_n_dims; i<n_dims; i++){
    inds_arr[i] = offsets[i-(n_dims-offsets_n_dims)];
  }
  //Now convert N-d position into 1-d index and do the copy
  indx = get_index(n_dims, inds_arr);
  if(indx==-1) return 1;
  std::copy(dat_in, dat_in + sz_in, data+indx);

  free(inds_arr);
  return 0;

}

bool my_array::populate_complex_slice(my_type * dat_in, size_t n_dims_in, size_t * offsets, size_t* sizes){
/** \brief Populate a slice of array from the input array.
*
* Extends populate_slice to fill an array slice from a larger array. As populate slice, assumes destination is a section of dimension m, with some finite offset in dimensions from m to n only. Assuming dat_in is a 1-d array in Fortran order (see get_element) this will read the proper subsection. Note this works fine for simple slices but costs more
 @param dat_in pointer to data 
 @param n_dims Dimensionality of input (must be less than dimension of array) 
 @param offsets Offsets in the other dimensions 
 @param sizes Sizes of input array 
 @return 0 (success), 1 else
 \todo Add testing of this
*/
  std::cerr<<"populate_complex_slice IS AN UNTESTED ROUTINE AT PRESENT !!!!!!!!!!!!!!!"<<'\n';
  if(n_dims_in >= n_dims) return 1;

  long indx;
  //Infer size of dat_in as product of sizes of all dims not sliced on
  size_t sz_in=1;
  size_t * inds_arr;
  inds_arr = (size_t *) calloc(n_dims,sizeof(size_t));
  for(size_t i=0; i<n_dims-n_dims_in; i++){
    sz_in *=dims[i];
  }
  //Location to copy to is 0 in all the matched dimensions and matches offsets in all the rest
  for(size_t i = n_dims-n_dims_in; i<n_dims; i++){
    inds_arr[i] = offsets[i-(n_dims-n_dims_in)];
  }
  //Now convert N-d position into 1-d index

  indx = get_index(n_dims, inds_arr);
  if(indx==-1) return 1;

  //Match up the sections and copy into place
  size_t input_n_dims= n_dims-n_dims_in, tot_sz=this->get_total_elements();
  std::vector<size_t> dest_inds_vec;
  size_t dest_ind, input_ind;
  for(size_t i=0; i< tot_sz; i++){
    //Destination index
    dest_ind = i + indx;
    dest_inds_vec = get_indices_from_offset(dest_ind);
    //Source index
    input_ind = dest_inds_vec[0];
    for(size_t j=1; j < input_n_dims; j++) input_ind += sizes[j-1]*dest_inds_vec[i];

    std::copy(dat_in + input_ind, dat_in+input_ind +1, data + dest_ind);
  }
  //Clumsy but hey

  return 0;

}

//Dummy to document file format
/**
    * @class dummy_file_format
    * The layout is:
    * sizeof(size_t) sizeof(my_type) io_verification_code Version string
    *Next_block n_dims dims[n_dims]
    *Next_block data
  *IMPORTANT: the VERSION specifier links output files to code. If the file output is changed, commit and clean build with a bumped major version number tag to correctly specify this
*/
bool my_array::write_to_file(std::fstream &file){
/** \brief Write array to file
*
*Writes array to file. Data is in a few blocks each starting with a number defining their end position in the file. \copydoc dummy_file_format
@param file Filestream to write to
@return 0 (success) 1 (error)
*/
  if(!file.is_open() || (this->data ==nullptr)) return 1;

  int size_sz = sizeof(size_t);
  int my_sz = sizeof(my_type);
  bool write_err = 0;
  
  //Next location is current position plus sizes we're about to write
  size_t next_location = (size_t) file.tellg() + 2*sizeof(int)+my_sz + GIT_VERSION_SIZE*sizeof(char);
  
  //Write Size of ints and data, IO verification constant and Code version.
  const char tmp_vers[GIT_VERSION_SIZE] = VERSION;
  file.write((char*) & size_sz, sizeof(int));
  file.write((char*) & my_sz, sizeof(int));
  file.write((char*) &io_verify, my_sz);
  file.write((char*) &tmp_vers, sizeof(char)*GIT_VERSION_SIZE);
 
  //Check file location matches what we expected
  if((size_t)file.tellg() != next_location) write_err=1;
  

  //Next prepare dimension info
  size_t total_size = get_total_elements();
  next_location += size_sz*(n_dims+2);
  //Write position of next section
  file.write((char*) & next_location, size_sz);
  //Now write the dims stuff
  size_t n_dims_tmp = n_dims;
  file.write((char*) & n_dims_tmp, size_sz);
  size_t dim_tmp;
  for(size_t i=0;i<n_dims;i++){
    dim_tmp = dims[i];
    file.write((char*) &dim_tmp, size_sz);
  }
  //Check file location matches what we expected
  if((size_t)file.tellg() != next_location) write_err=1;

  //Next we write the data
  next_location += my_sz*total_size + size_sz;
  file.write((char*) & next_location, size_sz);
  file.write((char *) data , my_sz*total_size);
  if((size_t)file.tellg() != next_location) write_err=1;

  //Report if any of our offsets didn't match the file
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);

  return 0;

}

bool my_array::write_section_to_file(std::fstream &file, std::vector<size_t> bounds){
/**\brief Write a subsection of array to file
*
*We write only the section delimited by bounds, which should have two elements for each dimension of this array. \copydoc dummy_file_format
@param file Filestream to write to
@param bounds Vector of indices delimiting subsection to write
@return 0 (success) 1 (error)
*/

  if(!file.is_open() || (this->data ==nullptr)) return 1;
  //Check we have enough bounds and they're all valid
  if(bounds.size() != n_dims*2) return 1;
  for(size_t i=0; i< n_dims; i++){
    if(bounds[2*i+1] > dims[i]) return 1;
  }
  
  int size_sz = sizeof(size_t);
  int my_sz = sizeof(my_type);
  bool write_err = 0;

//This section matches write_to_file exactly -------------------
  size_t next_location = (size_t) file.tellg() + 2*sizeof(int)+my_sz + GIT_VERSION_SIZE*sizeof(char);
  
  const char tmp_vers[GIT_VERSION_SIZE] = VERSION;
  file.write((char*) & size_sz, sizeof(int));
  file.write((char*) & my_sz, sizeof(int));
  file.write((char*) &io_verify, my_sz);
  file.write((char*) &tmp_vers, sizeof(char)*GIT_VERSION_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
//--------------------------------------------------------------

  //Write dims information, correcting the sizes to that of the subsection as we do
  size_t total_size = 1;
  next_location += size_sz*(n_dims+2);
  file.write((char*) & next_location, size_sz);
  
  size_t n_dims_tmp = n_dims;
  file.write((char*) & n_dims_tmp, size_sz);

  size_t dim_tmp;
  for(size_t i=0;i<n_dims;i++){
    dim_tmp = (bounds[i*2+1]-bounds[i*2]);
    file.write((char*) &dim_tmp, size_sz);
    total_size*=dim_tmp;
  }

  if((size_t)file.tellg() != next_location) write_err=1;

  //Now we write the data
  //We use get_element and nested loops for simplicity
  next_location += my_sz*total_size + size_sz;
  file.write((char*) & next_location, size_sz);

  my_type element;

  //Write each element to the file
  size_t * positions;
  positions = (size_t *) malloc(n_dims*sizeof(size_t));
  for(size_t dim = 0; dim < n_dims; dim ++) positions[dim] = bounds[dim*2];

  for(size_t el = 0; el < total_size; el++){
    element = get_element(n_dims, positions);
    file.write((char *) &element , sizeof(my_type));

    positions[0]++;//Increment lowest dim
    for(size_t dim = 0; dim < n_dims; dim ++){
      //Cascade increment up ensuring each dim stays in bounds
      if(positions[dim] >= bounds[dim*2 + 1]){
        positions[dim] = bounds[dim*2];
        if(dim < n_dims-1) positions[dim+1]++; //Only if not already the top dimension
      }else{
        //Stop looping when cascade ends
        break;
      }
    }
  }
  free(positions);

  if((size_t)file.tellg() != next_location) write_err=1;

  //Report if any of our offsets didn't match the file
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);

  return 0;

}

bool my_array::read_from_file(std::fstream &file){
/** \brief Read array from file
*
*Reads data from file. This array should have already been created in the correct shape, otherwise we return an error.
  \copydoc dummy_file_format
@param file Filestream to read from
@return 0 (success), 1 else
*/

  //Read the dimensions from file and check they match this array
  std::vector<size_t> dims_vec = read_dims_from_file(file);
  if(dims_vec.size() !=n_dims){
    my_error_print("Dimensions do not match, aborting read", mpi_info.rank);
    return 1;
  }
  size_t tot_els =1;
  for(size_t i=0;i<dims_vec.size();i++){
    if(dims_vec[i] !=dims[i]){
      my_error_print("Dimensions do not match, aborting read", mpi_info.rank);
      return 1;
    }
    tot_els *= dims[i];
  }

  //Read the next-block tag and then the data
  size_t next_block;
  file.read((char*) &next_block, sizeof(size_t));
  file.read((char *) data , sizeof(my_type)*tot_els);

  return 0;

}

std::vector<size_t> my_array::read_dims_from_file(std::fstream &file){
/** \brief Read dimensions from array file
*
*Reads dims from file into vector. Returns empty vector on read error \copydoc dummy_file_format
@param file Filestream to read from
@return Vector of dimensions
*/
  std::vector<size_t> dims_vec;
  char tmp_vers[GIT_VERSION_SIZE];
  my_type verf=0.0;

  if(!file.good()){
    my_error_print("File access error");
    return dims_vec;
  }

  size_t n_dims_in, dim_tmp, next_block;
  //Read the data size specifiers and check they match ours, otherwise error
  int size_sz=0, my_sz=0;
  file.read((char*) &size_sz, sizeof(int));
  file.read((char*) &my_sz, sizeof(int));

  if(size_sz !=sizeof(size_t)) my_error_print("size_t size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type)) my_error_print("my_type size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type) ||size_sz !=sizeof(size_t)) return dims_vec;
  
  //Read the verification constant and the Version code
  file.read((char*) &verf, sizeof(my_type));
  file.read((char*) &tmp_vers, sizeof(char)*GIT_VERSION_SIZE);

  if(verf != io_verify){
  //equality even though floats as should be identical
    //If we have a read error and the versions don't match, that might be the problem...
    my_error_print("File read error", mpi_info.rank);
    if(compare_as_version_string(tmp_vers, VERSION, true) != 0  && compare_as_version_string(tmp_vers, "IDL data write") != 0) my_error_print("Incompatible code version or file", mpi_info.rank);
    return dims_vec;
  }else{
    //Otherwise we report on verson mismatch according to considering major only
    if(compare_as_version_string(tmp_vers) != 0 && compare_as_version_string(tmp_vers, "IDL data write") != 0){
       my_print("WARNING: A different code version was used to write this data. Proceed with caution. Fields may not align correctly and data may differ", mpi_info.rank);
      std::string tmp = tmp_vers;
      my_print(tmp, mpi_info.rank);
      tmp = VERSION;
      my_print(tmp, mpi_info.rank);
    }
  }
  //Read the dimensions information
  file.read((char*) &next_block, sizeof(size_t));
  file.read((char*) &n_dims_in, sizeof(size_t));

  for(size_t i=0;i<n_dims_in;i++){
    file.read((char*) &dim_tmp, sizeof(size_t));
    dims_vec.push_back(dim_tmp);
  }

  return dims_vec;
}

/********Helpful functions working on entire array as a thing ****/

bool my_array::resize(size_t dim, size_t sz, bool verbose){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Note due to using 1-d memory layout both cases require copying all data and therefore briefly memory to store the old and new arrays. However shinking the last dimension does not necessarily require a copy.
@param dim Dimension to resize
@param sz New size for dimension
@param verbose Flag to print extra info
@return 0 (nothing to report) 1 else, including "new size matches, nothing done"
*/

  if(verbose) my_print("Attempting to resize", mpi_info.rank);

  if(sz == 0 || sz > MAX_SIZE) return 1;
  //size errors. We don't allow to resize any dimension to 0!!!
  if(dim > this->n_dims) return 1;
  //Shortcut if nothing to do
  if(sz == dims[dim]){
     if(verbose) my_print("Size matches", mpi_info.rank);
     return 1;
  }
  
  my_type * new_data;
  size_t part_sz = 1;

  if(dim == n_dims-1){
    //special case as we can shrink and maybe grow without copy
    for(size_t i=0; i<n_dims-1; ++i) part_sz*= dims[i];
    //Calculate product of all other dims

    new_data = (my_type *) realloc((void*) this->data, part_sz*sz*sizeof(my_type));
    if(!new_data){
      my_error_print("Failed to reallocate memory", mpi_info.rank);
      return 1;
      //failure. leave as was.
    }
    //zero the newly added elements
    long new_els = part_sz*(sz - dims[dim]);
    if(new_els > 0) memset((void*)(new_data + part_sz*dims[dim]), 0.0, new_els*sizeof(my_type));
  }
  else{
    //have to allocate a new block and copy across.
    //We do this in segments rather than element-wise.
    size_t els_to_copy = 1, n_segments = 1;

    for(size_t i=0; i<dim; ++i) els_to_copy *= dims[i];
    for(size_t i=dim+1; i< n_dims; ++i) n_segments *= dims[i];

    part_sz = els_to_copy*n_segments;
    new_data=(my_type*)calloc(part_sz*sz,sizeof(my_type));
    if(!new_data){
      my_error_print("Failed to allocate memory", mpi_info.rank);
      return 1;
      //failure. leave as was.
    }
    // Now we know dim is zeroth or a middle dimension. (1 for n_dims==3, 1 or 2 for n_dims==4.So we copy in chunks
    size_t chunk_sz = els_to_copy;
    (sz> dims[dim])? els_to_copy *= dims[dim] : els_to_copy *= sz;
    for(size_t i=0; i< n_segments; ++i) std::copy(data + i*chunk_sz*dims[dim],data + i*chunk_sz*dims[dim]+ els_to_copy, new_data + i*chunk_sz*sz);

    free(data);
  }
  //Set data to be the amended/new array and update dims
  data = new_data;
  dims[dim] = sz;

  if(verbose) my_print("New size of dim "+mk_str(dim)+  " is " + mk_str(dims[dim]), mpi_info.rank);

  return 0;

}

bool my_array::shift(size_t dim, long n_els){
/** \brief Shift array on dimension dim by n_els
*
* Because of individual getter/setter per dimensionality, we use the 1-d backing to do this.
* @param dim Dimension to shift, (0 to n_dims-1) 
@param n_els Number of elements to shift by. 
@return 0 (success) 1 else
\todo Fix special case
*/
/*
* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>
*/

  if(dim > this->n_dims) return 1;
  if(n_els ==0) return 0;
  //Correct n_els for being more than one full cycle and move < 0 to equivalent +ve
  long sign_n = (n_els<0? -1: 1);
  n_els = sign_n >0? (std::abs(n_els)%dims[dim]):dims[dim]-(std::abs(n_els)%dims[dim]) ;
  //By now n_els guaranteed >= 0
  
  my_type * new_data;
  size_t sub_sz = 1;

/*  if(dim == n_dims -1){
    size_t part_sz =0;
    //Special case else we'd be copying up to the whole array at once
    //We extract one chunk, hold it aside while we slot the rest into place and then put it back
    for(size_t i=0; i<dim-1; ++i) sub_sz*= dims[i];
    //Size of sub chunk which stays intact as we rotate

    size_t  n_segments = 1;
    size_t n_cyc=1;
    //Total size of copyable chunk
    size_t chunk_sz = sub_sz*dims[dim-1];
    new_data=(my_type*)malloc(chunk_sz*sizeof(my_type));
    size_t actual_shift = 0;

    //Total number of chunks
//    for(size_t i=dim; i< n_dims; ++i) n_segments *= dims[i];
    n_segments = dims[dim];
    n_cyc = n_segments | n_els;
    

    //Put one chunk aside safely
    size_t removed_chunk = 0;
    std::copy(data + removed_chunk*chunk_sz,data + (removed_chunk+1)*chunk_sz, new_data);
    size_t last_pos=0, dest_chunk =0;
    for(size_t i=0; i< n_segments-1; ++i){
      
      //Now move the chunk that replaces it
      dest_chunk = (last_pos-n_els >= dims[dim])? last_pos-n_els-dims[dim]: last_pos-n_els;
      dest_chunk += (dest_chunk < 0 ? dims[dim]:0);
      //copy the dest chunk to the last vacated position
      std::cout<<"working on "<<(*(data+last_pos*chunk_sz))<<'\n';
      for(size_t i=0; i<this->get_dims(0); i++){
        for(size_t j =0; j<this->get_dims(1); j++){
          std::cout<<this->get_element(i, j)<<" ";
        }
      std::cout<<'\n';
      }

      std::copy(data+dest_chunk*chunk_sz, data+(dest_chunk+1)*chunk_sz, data+last_pos*chunk_sz);
      last_pos=dest_chunk;
    }
    //Slot the spare chunk back in
    std::copy(new_data, new_data+chunk_sz, data+last_pos*chunk_sz);
    //Damn that only works for relatively prime size and swap
    //This uses as many element copies as the other case version but less memory
    for(size_t i=0; i< n_segments-1; ++i){
      //identify which chunk lives here
      //Now move the chunk that replaces it
      dest_chunk = (last_pos-n_els >= dims[dim])? last_pos-n_els-dims[dim]: last_pos-n_els;
      dest_chunk += (dest_chunk < 0 ? dims[dim]:0);
      for(size_t i=0; i<this->get_dims(0); i++){
        for(size_t j =0; j<this->get_dims(1); j++){
          std::cout<<this->get_element(i, j)<<" ";
        }
      std::cout<<'\n';
      }

      std::copy(data+dest_chunk*chunk_sz, data+(dest_chunk+1)*chunk_sz, data+last_pos*chunk_sz);
      last_pos=dest_chunk;
    }
    
  }else*/
  {
    //allocate a new block and copy across. chunking so we don't need to use element getters
    for(size_t i=0; i<dim; ++i) sub_sz*= dims[i];
    //Size of sub chunk which stays intact as we rotate

    size_t  n_segments = 1;
    //Total size of copyable, rotateable chunk
    size_t chunk_sz = sub_sz*dims[dim];
    new_data=(my_type*)malloc(chunk_sz*sizeof(my_type));
    long actual_shift = 0;

    //Total number of chunks
    for(size_t i=dim+1; i< n_dims; ++i) n_segments *= dims[i];
    
    //Now we rotate each chunk
    for(size_t i=0; i< n_segments; ++i){
      //We extract the first chunk to our spare memory
      std::copy(data + i*chunk_sz,data + (i+1)*chunk_sz, new_data);
//      std::cout<<*(data+i*chunk_sz+2)<<" "<<*(new_data+2)<<'\n';
          //And rotate it as we put it back. This is two copies, but damn is it easier
      for(size_t j=0; j< dims[dim]; j++){
        actual_shift = (j+n_els >= dims[dim])? j+n_els-dims[dim]: j+n_els;
        actual_shift += (actual_shift < 0 ? dims[dim]:0);
        //std::cout<<actual_shift<<" "<<j+n_els<<" "<<j+n_els-dims[dim]<<'\n';
        std::copy(new_data+j*sub_sz, new_data+(j+1)*sub_sz, data + i*chunk_sz+actual_shift*sub_sz);
      }
    }
  }
  if(new_data) free(new_data);

  return 0;

}

my_type my_array::minval(size_t offset){
/** \brief Find minimum value of data
*
*Finds the minimum of the array after offset, using linear search through contiguous memory. Offset should be obtained using one of the get_index functions.
@param offset 1-D array offset to start search from
@return The minimum value of data
*/
  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::min();

  return *(std::min_element(data+offset, data+total_size));

}
my_type my_array::minval(std::vector<size_t> &ind, size_t offset){
/** \brief Find minimum value of data
*
*Finds the minimum of the array after offset, using linear search through contiguous memory. Offset should be obtained using one of the get_index functions.
@param offset 1-D array offset to start search from
@param[out] ind The indices where min value is located
@return The minimum value of data
*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::min();

  auto it = std::min_element(data+offset, data+total_size);
  ind = get_indices_from_offset(it - data);
  return *(it);
  
}

my_type my_array::maxval(size_t offset){
/** \brief Find maximum value of data
*
*Finds the maximum of the array after offset, using linear search through contiguous memory. Offset should be obtained using one of the get_index functions.
@param offset 1-D array offset to start search from
@return The maximum value of data
*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  return *(std::max_element(data+offset, data+total_size));

}
my_type my_array::maxval(std::vector<size_t> &ind, size_t offset){
/** \brief Find maximum value of data
*
*Finds the maximum of the array after offset, using linear search through contiguous memory. Offset should be obtained using one of the get_index functions.
@param offset 1-D array offset to start search from
@param[out] ind The indices where max value is located
@return The maximum value of data
*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  auto it = std::max_element(data+offset, data+total_size);
  ind = get_indices_from_offset(it - data);
  return *(it);

}

my_type my_array::partial_maxval(std::vector<std::pair<size_t, size_t> > ranges, std::vector<size_t> &ind){
/** \brief Maximum value over part of range
*
* WARNING: this is slow. Perhaps very slow. I'm using routines I have to knock it up quickly. Beware!!!
@param ranges The indices to consider between on each dimension
@param[out] ind The indices where max was located
@return The maximum value over given ranges
*/

  if(ranges.size() != this->n_dims) return std::numeric_limits<my_type>::max();

  my_array tmp = *(this);
  //Work on a copy
  for(size_t i=0; i< ranges.size(); i++){
    tmp.shift(i, -ranges[i].first);
    tmp.resize(i, (ranges[i].second-ranges[i].first));
  }
  my_type val = tmp.maxval(ind);
  for(size_t i=0; i< ranges.size(); i++) ind[i] +=ranges[i].first;
  return val;

}

void my_array::smooth_1d(int n_pts){
/** \brief Smooth a 1-d data_array
*
*Smooths the data backing array. This does strange things for non-1d data at the ends.
@param n_pts Smoothing width
*/
  if(n_dims !=1) return;
  inplace_boxcar_smooth(data, (int) get_total_elements(), n_pts);

}
