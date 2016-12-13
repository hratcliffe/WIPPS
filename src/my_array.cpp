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

extern const mpi_info_struc mpi_info;

my_array::my_array(){
/** Default constructor*/
  construct();
}

void my_array::construct(){
/** \brief Shared contructor code
*
*Sets default values of things
*/
  n_dims = 0;
  data=nullptr;
  dims=nullptr;

}

void my_array::alloc_all(const size_t n_dims, const size_t * const dims){
/** \brief Takes care of memory allocation
*
* Allocate dims array and data. If any size is zero, any exceeds MAX_SIZE, or overall exceeds MAX_SIZE_TOT, print error and stop. NB inputs will not be modified, will be copied
*/

  size_t tot_dims = 1;
  bool too_large=false;
  for(size_t i=0; i< n_dims; i++){
    tot_dims*=dims[i];
    if(dims[i] > MAX_SIZE) too_large = true;
  }
  //tot_dims now 0 if any dim is 0
  if(tot_dims==0){
    my_print("Array cannot have 0 dim", mpi_info.rank);
    return;
  }
  if(too_large || tot_dims > MAX_SIZE_TOT){
    if(too_large) my_print("Array size exceeds max per-dimension size of "+mk_str(MAX_SIZE), mpi_info.rank);
    if(tot_dims> MAX_SIZE_TOT) my_print("Array size exceeds max overall size of "+mk_str(MAX_SIZE_TOT), mpi_info.rank);
    return;
  }
  
  this->n_dims =n_dims;
  this->dims = (size_t*)malloc(n_dims*sizeof(size_t));
  if(dims){
  //If allocation succeeds, we can copy in dims data and allocate data array, else stop
    std::copy(dims, dims+n_dims, this->dims);
    data=(my_type*)calloc(tot_dims,sizeof(my_type));
  }

}
my_array::my_array(size_t nx, size_t ny, size_t nz, size_t nt){
/** \brief 1 to 4 d rectangular array helper
*
*Sets up a n-d rectangular array for n = 1 to 4. Helper avoids user having to construct size_t array of dims
*/
  construct();
  
  size_t * dims_in;
  size_t n_dims_in, n_dims_act = 4;
  n_dims_in = 4;
  dims_in = (size_t*)malloc(n_dims_in*sizeof(size_t));
  dims_in[0] = nx;
  dims_in[1] = ny;
  dims_in[2] = nz;
  dims_in[3] = nt;
  
  for(size_t i=4; i>0; i--){
    if(dims_in[i-1] == 0) n_dims_act --;
    else break;
  }
  //Check which dims are zero
  alloc_all(n_dims_act, dims_in);
  
  free(dims_in);
  
}

my_array::my_array(size_t n_dims, size_t * dims ){
/** \brief Arbitrary dim rectangular array
*
*Sets up a n_dims array of dimensions dims and allocates data arrays.
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

my_array & my_array::operator=(const my_array& src){
/** \brief Copy assignment
*
*Sets this equal to a copy of source
*/


  if(this->data) free(data);
  if(this->dims) free(dims);
  //Clean up in case this was already a thing
  my_array::construct();

  if(!src.is_good()) return *this;
  //No copy if src is zero
  
  alloc_all(src.n_dims, src.dims);
  size_t tot_els = this->get_total_elements();
  if(this->data && src.data) std::copy(src.data, src.data + tot_els, this->data);

  return *this;

}

my_array::my_array(my_array && src){
/** \brief Move constructor
*
*Move src to a new instance i.e. copy fields but don't move memory. Src becomes empty afterwards
*/

  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims
  
  this->n_dims = src.n_dims;
  src.n_dims = 0;
  this->dims = src.dims;
  src.dims = nullptr;
  this->data = src.data;
  src.data = nullptr;
  //Steal memory from src
}

my_array::my_array(const my_array &src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data
*/

  construct();
  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims
  
  alloc_all(src.n_dims, src.dims);
  size_t tot_els = this->get_total_elements();
  if(this->data && src.is_good()) std::copy(src.data, src.data + tot_els, this->data);

}

my_type * my_array::disown_data(){
/** \brief Disown and return data pointer
*
*Surrenders ownership of memory pointed to by data, nullifies dimensions and returns pointer. NB if this pointer is not kept an manually freed, memory will leak
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
* Duplicates sizes of src but do not copy data
*/
  construct();

  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims

  //Delete any previous allocation
  if(data) free(data);
  if(dims) free(dims);

  alloc_all(src.n_dims, src.dims);
  
}

long my_array::get_index(size_t n_dims, size_t * inds_in)const{
/** Passing integer array and loop so going to be slower
*/
  if(n_dims != this-> n_dims) return -1;
  for(size_t i=0; i< n_dims; i++){
   // std::cout<<i<<' '<<inds_in[i]<<' ';
      if(inds_in[i] >=dims[i]) return -1;
  }

  long ind = inds_in[0];
  size_t subdims = dims[0];
  for(size_t i=1; i< n_dims; i++){
    ind += subdims*inds_in[i];
    subdims *= dims[i];
  }
  return ind;
  
}
long my_array::get_index(size_t nx)const{
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args

* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>
Etc */

  if(n_dims != 1){
#ifdef DEBUG_DIMS
    my_print("Wrong array dimension, attempting 1 with "+mk_str(n_dims), mpi_info.rank);
#endif
    return -1;

  }
  if((nx < dims[0])){
    return (long) nx;
  }else{
    return -1;
  }

}
long my_array::get_index(size_t nx, size_t ny)const{
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args

* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is <-row->
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<---------slice------->
Etc
*/

  if(n_dims != 2){
#ifdef DEBUG_DIMS
    my_print("Wrong array dimension, attempting 2 with "+mk_str(n_dims), mpi_info.rank);
#endif
    return -1;

  }
  if((nx < dims[0]) && (ny<dims[1])){
    return (long) (ny*dims[0] + nx);
  }else{
    return -1;
  }
}
long my_array::get_index(size_t nx, size_t ny, size_t nz)const{
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args
*/

  if(n_dims != 3){
#ifdef DEBUG_DIMS
    my_print("Wrong array dimension, attempting 3 with "+mk_str(n_dims), mpi_info.rank);
#endif
    return -1;
  }

  if((nx < dims[0]) && (ny<dims[1]) && ((nz<dims[2]))){
    return (nz*dims[1]+ ny)*dims[0] + nx;
  }else{
    return -1;
  }
}
long my_array::get_index(size_t nx, size_t ny, size_t nz, size_t nt)const{
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args
*/

  if(n_dims != 4){
#ifdef DEBUG_DIMS
    my_print("Wrong array dimension, attempting 4 with "+mk_str(n_dims), mpi_info.rank);
#endif
    return -1;
  }
  if((nx < dims[0]) && (ny<dims[1])&& ((nz<dims[2]))&& ((nt<dims[3]))){
    return ((nt*dims[2] +nz)*dims[1]+ ny)*dims[0] + nx;
  }else{
    return -1;
  }
}

std::vector<size_t> my_array::get_index_from_offset(size_t offset)const{

  std::vector<size_t> pos;
  pos.resize(n_dims);
  size_t tmp_offset = offset;
  
  size_t stride=1;
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

size_t my_array::get_dims()const{
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

my_type my_array::get_element(size_t nx)const{
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(nx);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t nx, size_t ny)const{
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz)const{
/** Return element at nx, ny, nz. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny, nz);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz, size_t nt)const{
/** Return element at nx, ny, nz, nt. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny, nz, nt);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t n_dims, size_t* dims)const{
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(n_dims, dims);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element_from_index(size_t ind)const{
  if(ind < get_total_elements()){
    return *(data+ind);
  }else{
    return 0.0;
  }
}

size_t my_array::get_total_elements()const{
/** Return total size of array */
  size_t tot_els=1;

  for(size_t i=0; i<n_dims;++i) tot_els *=dims[i];
  return tot_els;

}

bool my_array::set_element(size_t nx, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, @return 1 if out of range, wrong number of args, 0 else.
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
/** \brief Sets array element
*
*Sets elements at nx, ny, @return 1 if out of range, wrong number of args, 0 else.
*/

  long index = get_index(nx, ny);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}
bool my_array::set_element(size_t nx, size_t ny, size_t nz, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, @return 1 if out of range, wrong number of args, 0 else.
*/

  long index = get_index(nx, ny, nz);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}
bool my_array::set_element(size_t nx, size_t ny, size_t nz, size_t nt, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, nt, @return 1 if out of range, wrong number of args, 0 else.
*/

  long index = get_index(nx, ny, nz, nt);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}
bool my_array::set_element(size_t n_dims_in, size_t * dims_in, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, nt, @return 1 if out of range, wrong number of args, 0 else.
*/

  long index = get_index(n_dims_in, dims_in);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  

}

bool my_array::populate_data(my_type * dat_in, size_t n_tot, bool convert){
/** \brief Fill array
*
*Populates data with the total number of elements specified, as long as n_tot is less than the product of dims. Parameter needed so we can't overflow dat_int. Assumes same row-column order etc etc @return 0 (sucess) 1 (error). If convert is set we convert to the other of double and float in the copy
*/

  size_t tot_els = get_total_elements();
//  if(n_tot > tot_els) return 1;
  if(n_tot < tot_els) tot_els = n_tot;
  //Use min of n_tot, tot_els


  if(!convert) std::copy(dat_in, dat_in+tot_els, this->data);
  else std::copy( (other_type *) dat_in, (other_type *) dat_in+tot_els, this->data);
  return 0;

}

bool my_array::populate_slice(my_type * dat_in, size_t n_dims_in, size_t * offsets){
/** \brief Populate a slice of the array from the input array.
*
* Extends populate_row to fill an array slice, that is a section of dimension m, with some finite offset in dimensions from m to n only. E.g. to fill a row of a 3-d array call with n_dims_in = 3-1=2 and offsets={column, plane}
 @param dat_in pointer to data @param n_dims Dimensionality of input (must be less than dimension of array) @param offsets Offsets in the other dimensions
*/

  if(n_dims_in >= n_dims) return 1;

  long indx;
  size_t sz_in=1;
  size_t * inds_arr;
  inds_arr = (size_t *) malloc(n_dims*sizeof(size_t));
  for(size_t i=0; i<n_dims-n_dims_in; i++){
    inds_arr[i] = 0;
    sz_in *=dims[i];
  }
  for(size_t i = n_dims-n_dims_in; i<n_dims; i++){
    inds_arr[i] = offsets[i-(n_dims-n_dims_in)];
  }
  indx = get_index(n_dims, inds_arr);
  if(indx==-1) return 1;

  std::copy(dat_in, dat_in + sz_in, data+indx);

  free(inds_arr);
  return 0;

}

bool my_array::populate_complex_slice(my_type * dat_in, size_t n_dims_in, size_t * offsets, size_t* sizes){
/** \brief Populate a slice of the array from the input array.
*
* Extends populate_slice to fill an array slice from a larger array. As populate slice, assumes destination is a section of dimension m, with some finite offset in dimensions from m to n only. Assuming dat_in is a 1-d array in Fortran order (see get_element) this will read the proper subsection. Note this works fine for simple slices but costs more
 @param dat_in pointer to data @param n_dims Dimensionality of input (must be less than dimension of array) @param offsets Offsets in the other dimensions @param sizes Sizes of input array
*/
  std::cout<<"populate_complex_slice IS AN UNTESTED ROUTINE!!!!!!!!!!!!!!!"<<'\n';
  if(n_dims_in >= n_dims) return 1;

  long indx;
  size_t sz_in=1;
  size_t * inds_arr;
  inds_arr = (size_t *) malloc(n_dims*sizeof(size_t));
  for(size_t i=0; i<n_dims-n_dims_in; i++){
    inds_arr[i] = 0;
    sz_in *=dims[i];
  }
  for(size_t i = n_dims-n_dims_in; i<n_dims; i++){
    inds_arr[i] = offsets[i-(n_dims-n_dims_in)];
  }
  indx = get_index(n_dims, inds_arr);

  if(indx==-1) return 1;

  size_t input_n_dims= n_dims-n_dims_in, tot_sz=this->get_total_elements();
  std::vector<size_t> dest_inds_vec;
  size_t dest_ind, input_ind;
  for(size_t i=0; i< tot_sz; i++){
    //Destination index
    dest_ind = i + indx;
    dest_inds_vec = get_index_from_offset(dest_ind);
    //Source index
    input_ind = dest_inds_vec[0];
    for(size_t j=1; j < input_n_dims; j++) input_ind += sizes[j-1]*dest_inds_vec[i];

    std::copy(dat_in + input_ind, dat_in+input_ind +1, data + dest_ind);
  }
  //Clumsy but hey

  return 0;

}

bool my_array::write_to_file(std::fstream &file){
/** \brief Write array to file
*
*Writes array to file. Data is in a few blocks each starting with a number defining their end position in the file. The layout is:
* Next_block sizeof(size_t) sizeof(my_type) io_verification_code Version string
*Next_block n_dims dims[n_dims]
*Next_block data
  *IMPORTANT: this VERSION specifier links output files to code. If modifying output or order commit and clean build before using. @return 0 (sucess) 1 (error)
*/
  if(!file.is_open() || (this->data ==nullptr)) return 1;

  int size_sz = sizeof(size_t);
  int my_sz = sizeof(my_type);
  bool write_err = 0;
  
  size_t next_location = (size_t) file.tellg() + 2*sizeof(int)+my_sz + 15*sizeof(char);
  
  const char tmp_vers[15] = VERSION;
  file.write((char*) & size_sz, sizeof(int));
  file.write((char*) & my_sz, sizeof(int));
  file.write((char*) &io_verify, my_sz);
  file.write((char*) &tmp_vers, sizeof(char)*15);
  //Sizeof ints and data, IO verification constant and Code version...
  std::cout<<"here"<<'\n';
std::cout<<size_sz<<' '<<my_sz<<' '<<io_verify<<'\n';
  if((size_t)file.tellg() != next_location) write_err=1;
  
  size_t total_size = get_total_elements();
  //dimension info
  next_location += size_sz*(n_dims+2);
  file.write((char*) & next_location, size_sz);
  //Position of next section
  
  size_t n_dims_tmp = n_dims;
  file.write((char*) & n_dims_tmp, size_sz);
  size_t dim_tmp;
  for(size_t i=0;i<n_dims;i++){
    dim_tmp = dims[i];
    file.write((char*) &dim_tmp, size_sz);
  }

  if((size_t)file.tellg() != next_location) write_err=1;

  next_location += my_sz*total_size + size_sz;
  file.write((char*) & next_location, size_sz);
  //Position of next section

  file.write((char *) data , my_sz*total_size);

  if((size_t)file.tellg() != next_location) write_err=1;

  if(write_err) my_print("Error writing offset positions", mpi_info.rank);

  return 0;

}

bool my_array::write_section_to_file(std::fstream &file, std::vector<size_t> bounds){
/**\brief Write a subsection of array to file
*
*See write_to_file for format

  * We use lazy method of get_element for each element, less prone to offset errors and memory is already much faster than disk
  *
  *IMPORTANT: this VERSION specifier links output files to code. If modifying output or order commit and clean build before using. @return 0 (sucess) 1 (error) \todo Probably need arb dims version...
*/

  if(!file.is_open() || (this->data ==nullptr)) return 1;
  if(bounds.size() != n_dims*2) return 1;
  for(size_t i=0; i< n_dims; i++){
    if(bounds[2*i+1] > dims[i]) return 1;
  }
  
  int size_sz = sizeof(size_t);
  int my_sz = sizeof(my_type);
  bool write_err = 0;
  
  size_t next_location = (size_t) file.tellg() + 2*sizeof(int)+my_sz + 15*sizeof(char);
  
  const char tmp_vers[15] = VERSION;
  file.write((char*) & size_sz, sizeof(int));
  file.write((char*) & my_sz, sizeof(int));
  file.write((char*) &io_verify, my_sz);
  file.write((char*) &tmp_vers, sizeof(char)*15);
  //Sizeof ints and data, IO verification constant and Code version...

  if((size_t)file.tellg() != next_location) write_err=1;

  size_t total_size = 1;
  //dimension info
  next_location += size_sz*(n_dims+2);
  file.write((char*) & next_location, size_sz);
  //Position of next section
  
  size_t n_dims_tmp = n_dims;
  file.write((char*) & n_dims_tmp, size_sz);

  size_t dim_tmp;
  for(size_t i=0;i<n_dims;i++){
    dim_tmp = (bounds[i*2+1]-bounds[i*2]);
    file.write((char*) &dim_tmp, size_sz);
    total_size*=dim_tmp;
  }

  if((size_t)file.tellg() != next_location) write_err=1;

  next_location += my_sz*total_size + size_sz;
  file.write((char*) & next_location, size_sz);
  //Position of next section

  my_type element;
  if(n_dims ==1){
    for(size_t i= bounds[0]; i< bounds[1]; i++){
      element = get_element(i);
      file.write((char *) &element, sizeof(my_type));
    }
  }else if(n_dims ==2){
    for(size_t j= bounds[2]; j< bounds[3]; j++){
      for(size_t i= bounds[0]; i< bounds[1]; i++){
        element = get_element(i, j);
        file.write((char *)  &element, sizeof(my_type));
      }
    }
  }else if(n_dims ==3){
    for(size_t k= bounds[4]; k< bounds[5]; k++){
      for(size_t j= bounds[2]; j< bounds[3]; j++){
        for(size_t i= bounds[0]; i< bounds[1]; i++){
          element  = get_element(i, j, k);
          file.write((char *) &element, sizeof(my_type));
        }
      }
    }
  }else if(n_dims ==4){
    for(size_t l= bounds[6]; l< bounds[7]; l++){
      for(size_t k= bounds[4]; k< bounds[5]; k++){
        for(size_t j= bounds[2]; j< bounds[3]; j++){
          for(size_t i= bounds[0]; i< bounds[1]; i++){
            element = get_element(i, j, k, l);
            file.write((char *) &element , sizeof(my_type));
          }
        }
      }
    }
  }
  
  if((size_t)file.tellg() != next_location) write_err=1;

  if(write_err) my_print("Error writing offset positions", mpi_info.rank);

  return 0;

}

bool my_array::read_from_file(std::fstream &file, bool no_version_check){
/** \brief Read array dump file
*
*Reads dims and data from file. Requires the dimensions of array to be already setup and will check for consistency with those in file

*/
/** sizeof(size_t) sizeof(my_type) io_verification_code Version string
*Next_block n_dims dims[n_dims]
*Next_block data
*/

  std::vector<size_t> dims_vec = read_dims_from_file(file, no_version_check);

  size_t next_block;
  
  if(dims_vec.size() !=n_dims){
    my_print("Dimensions do not match, aborting read", mpi_info.rank);
    return 1;
  }
  size_t tot_els =1;
  for(size_t i=0;i<dims_vec.size();i++){
    if(dims_vec[i] !=dims[i]){
      my_print("Dimensions do not match, aborting read", mpi_info.rank);
      return 1;
    }
    tot_els *= dims[i];
    //Check dims
  }
  file.read((char*) &next_block, sizeof(size_t));

  file.read((char *) data , sizeof(my_type)*tot_els);

  return 0;

}

std::vector<size_t> my_array::read_dims_from_file(std::fstream &file, bool no_version_check){
/** \brief Read array dump file for dimension
*
*Reads dims into vector. Returns empty vector on read error

*/
/** sizeof(size_t) sizeof(my_type) io_verification_code Version string
*Next_block n_dims dims[n_dims]
*Next_block data
*/
  std::vector<size_t> dims_vec;

  char tmp_vers[15];
  my_type verf=0.0;

  if(!file.good()){
    my_print("File access error");
    return dims_vec;
  }

  size_t n_dims_in, dim_tmp, next_block;
  int size_sz=0, my_sz=0;
  file.read((char*) &size_sz, sizeof(int));
  file.read((char*) &my_sz, sizeof(int));

  if(size_sz !=sizeof(size_t)) my_print("size_t size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type)) my_print("my_type size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type) ||size_sz !=sizeof(size_t)) return dims_vec;
  
  file.read((char*) &verf, sizeof(my_type));
  file.read((char*) &tmp_vers, sizeof(char)*15);

  if(verf != io_verify){
  //equality even though floats as should be identical
    my_print("File read error", mpi_info.rank);
    if(strcmp(tmp_vers, VERSION) !=0 && strcmp(tmp_vers, "IDL data write")!= 0 ) my_print("Incompatible code versions", mpi_info.rank);
    return dims_vec;
  }else{
    if(!no_version_check && strcmp(tmp_vers, VERSION) !=0){
       my_print("WARNING: A different code version was used to write this data. Proceed with caution. Fields may not align correctly.", mpi_info.rank);
      std::string tmp = tmp_vers;
      my_print(tmp, mpi_info.rank);
      tmp = VERSION;
      my_print(tmp, mpi_info.rank);
    }
  }

  file.read((char*) &next_block, sizeof(size_t));
  file.read((char*) &n_dims_in, sizeof(size_t));

  for(size_t i=0;i<n_dims_in;i++){
    file.read((char*) &dim_tmp, sizeof(size_t));
    dims_vec.push_back(dim_tmp);
  }

  return dims_vec;
}

bool my_array::resize(size_t dim, size_t sz, bool verbose){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Note due to using 1-d memory layout both cases require copying all data and therefore briefly memory to store the old and new arrays. However shinking the last dimension does not necessarily require a copy.
*/

  if(verbose) my_print("Attempting to resize", mpi_info.rank);

  if(sz == 0 || sz > MAX_SIZE)return 1;
  //size errors. Don't allow to resize to 0!!!
  if(dim > this->n_dims) return 1;
  if(sz == dims[dim]){
     if(verbose) my_print("Size matches", mpi_info.rank);
     return 1;
  }
  
  my_type * new_data;
  size_t part_sz = 1;

  if(dim == n_dims-1){
    //special case as we can shrink and maybe grow without copy
    for(size_t i=0; i<n_dims-1; ++i) part_sz*= dims[i];
    //product of all other dims

    new_data = (my_type *) realloc((void*) this->data, part_sz*sz*sizeof(my_type));
    if(!new_data){
      my_print("Failed to reallocate memory", mpi_info.rank);
      return 1;
      //failure. leave as was.
    }

    long new_els = part_sz*(sz - dims[dim]);
    
    if(new_els > 0) memset((void*)(new_data + part_sz*dims[dim]), 0.0, new_els*sizeof(my_type));
    //zero new elements
  }
  else{
    //have to allocate a new block and copy across.
    size_t els_to_copy = 1, n_segments = 1;

    for(size_t i=0; i<dim; ++i) els_to_copy *= dims[i];
    for(size_t i=dim+1; i< n_dims; ++i) n_segments *= dims[i];

//    for(size_t i=0; i<dim; ++i) part_sz*= dims[i];
    part_sz = els_to_copy*n_segments;
    new_data=(my_type*)calloc(part_sz*sz,sizeof(my_type));

    // Now we know dim is zeroth or a middle dimension. (1 for n_dims==3, 1 or 2 for n_dims==4.So we copy in chunks
    size_t chunk_sz = els_to_copy;
    (sz> dims[dim])? els_to_copy *= dims[dim] : els_to_copy *= sz;
    for(size_t i=0; i< n_segments; ++i) std::copy(data + i*chunk_sz*dims[dim],data + i*chunk_sz*dims[dim]+ els_to_copy, new_data + i*chunk_sz*sz);

    free(data);
  }

  data = new_data;
  dims[dim] = sz;

  if(verbose) my_print("New size of dim "+mk_str(dim)+  " is " + mk_str(dims[dim]), mpi_info.rank);

  return 0;

}

bool my_array::shift(size_t dim, long n_els){
/** \brief Shift array on dimension dim by n_els
*
* Because of individual getter/setter per dimensionality, we use the 1-d backing to do this.
* @param dim Dimension to shift, (0 to n_dims-1) @param n_els Number of elements to shift by. \todo Fix special case
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
/** Find minimum value of data, allows linear search through contiguous memory*/
  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::min();

  return *(std::min_element(data+offset, data+total_size));

}
my_type my_array::minval(std::vector<size_t> &ind, size_t offset){
/** Find minimum value of data, allows linear search through contiguous memory*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::min();

  auto it = std::min_element(data+offset, data+total_size);
  ind = get_index_from_offset(it - data);
  return *(it);
  

}

my_type my_array::maxval(size_t offset){
/** Find maximum value of data, allows linear search through contiguous memory*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  return *(std::max_element(data+offset, data+total_size));
}
my_type my_array::maxval(std::vector<size_t> &ind, size_t offset){
/** Find maximum value of data, allows linear search through contiguous memory*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  auto it = std::max_element(data+offset, data+total_size);
  ind = get_index_from_offset(it - data);
//  ind.push_back(it - data);

  return *(it);

}

my_type my_array::partial_maxval(std::vector<std::pair<size_t, size_t> > ranges, std::vector<size_t> &ind){
/** \brief Maximum value over part of range
*
* WARNING: this is slow. Perhaps very slow. I'm using routines I have to knock it up quickly. Beware!!!
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
my_type my_array::avval(){
/** Find average value of data */
  my_type av = 0.0;
  size_t total_size=get_total_elements();

  av=std::accumulate(data, data+total_size, 0.0);
  return av/total_size;
}

void my_array::smooth_1d(int n_pts){
  if(n_dims !=1) return;
  inplace_boxcar_smooth(data, (int) get_total_elements(), n_pts);

}
