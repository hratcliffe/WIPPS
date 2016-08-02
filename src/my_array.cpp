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
#include <algorithm>
#include <complex.h>
#include <limits.h>
#include <fftw3.h>
#include <cmath>
#include "support.h"
#include "my_array.h"

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
  defined = false;
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
  
  if(data && dims){
    defined = true;
  }
}

my_array::my_array(size_t n_dims, size_t * dims ){
/** \brief arbitrary dim rectangular array
*
*Sets up a n-d rectangular array and allocates data arrays.
*/
  construct();
  alloc_all(n_dims, dims);
  
  if(data && dims){
    defined = true;
  }
}

my_array::~my_array(){
/** Clean up explicit allocations
*
*
*/
  if(data) free(data);
  data=nullptr;
  if(dims) free(dims);
  dims=nullptr;

}

my_array & my_array::operator=(const my_array& src){
  
  if(this->data) free(data);
  this->construct();
  if(!src.dims || src.n_dims==0) return *this;
  //No copy if src is zero
  
  alloc_all(src.n_dims, src.dims);
  if(data && dims){
    defined = true;
  }
  size_t tot_els = this->get_total_elements();
  if(this->data && src.data) std::copy(src.data, src.data + tot_els, this->data);

  return *this;

}

my_array::my_array(const my_array &src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data \todo move constructor
*/
  construct();
  if(!src.dims || src.n_dims==0) return;
  //Stop if src has no dims
  
  alloc_all(src.n_dims, src.dims);
  if(data && dims){
    defined = true;
  }
  size_t tot_els = this->get_total_elements();
  if(this->data && src.is_good()) std::copy(src.data, src.data + tot_els, this->data);

}

long my_array::get_index(size_t n_dims, size_t * inds_in){
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
long my_array::get_index(size_t nx){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args

* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>
Etc \todo We may get speedup from removing checks. If so, wrap them in a debug IFDEF for fiddling vs running
*/

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
long my_array::get_index(size_t nx, size_t ny){
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
long my_array::get_index(size_t nx, size_t ny, size_t nz){
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
long my_array::get_index(size_t nx, size_t ny, size_t nz, size_t nt){
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

std::vector<size_t> my_array::get_index_from_offset(size_t offset){

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

size_t my_array::get_dims(){
  return n_dims;
}
size_t my_array::get_dims(size_t dim){
/** \brief Return size of dimension dim
*
*/
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

size_t my_array::get_length(size_t dim){
/** \brief Get size of dimension dim
*
**/
  return get_dims(dim);
}

my_type my_array::get_element(size_t nx){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(nx);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t nx, size_t ny){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz){
/** Return element at nx, ny, nz. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny, nz);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t nx, size_t ny, size_t nz, size_t nt){
/** Return element at nx, ny, nz, nt. Out of range etc will return 0.0*/
  long ind = get_index(nx, ny, nz, nt);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(size_t n_dims, size_t* dims){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  long ind = get_index(n_dims, dims);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}

size_t my_array::get_total_elements(){
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
bool my_array::set_element(size_t n_dims, size_t * dim, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, nt, @return 1 if out of range, wrong number of args, 0 else.
*/

  long index = get_index(n_dims, dims);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  

}

bool my_array::populate_data(my_type * dat_in, size_t n_tot){
/** \brief Fill array
*
*Populates data with the total number of elements specified, as long as n_tot is less than the product of dims. Parameter needed so we can't overflow dat_int. Assumes same row-column order etc etc @return 0 (sucess) 1 (error)
*/

  size_t tot_els = get_total_elements();
//  if(n_tot > tot_els) return 1;
  if(n_tot < tot_els) tot_els = n_tot;
  //Use min of n_tot, tot_els

  void * tmp = (void *) this->data;
  if(!tmp) return 1;

  memcpy (tmp , dat_in, tot_els*sizeof(my_type));

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
    for(int j=1; j < input_n_dims; j++) input_ind += sizes[j-1]*dest_inds_vec[i];

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

  if(file.tellg() != next_location) write_err=1;
  
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

  if(file.tellg() != next_location) write_err=1;

  next_location += my_sz*total_size + size_sz;
  file.write((char*) & next_location, size_sz);
  //Position of next section

  file.write((char *) data , my_sz*total_size);

  if(file.tellg() != next_location) write_err=1;

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

  if(file.tellg() != next_location) write_err=1;

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

  if(file.tellg() != next_location) write_err=1;

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
  
  if(file.tellg() != next_location) write_err=1;

  if(write_err) my_print("Error writing offset positions", mpi_info.rank);

  return 0;

}

bool my_array::read_from_file(std::fstream &file, bool no_version_check){
/** \brief Read array dump file
*
*Reads block id, data and axes from a file. Requires the dimensions of array to be already setup and will check for consistency with those in file

*/
/** sizeof(size_t) sizeof(my_type) io_verification_code Version string
*Next_block n_dims dims[n_dims]
*Next_block data
*/
  char tmp_vers[15];
  my_type verf=0.0;

  size_t n_dims_in, dim_tmp, next_block;
  int size_sz=0, my_sz=0;
  file.read((char*) &size_sz, sizeof(int));
  file.read((char*) &my_sz, sizeof(int));

  if(size_sz !=sizeof(size_t)) my_print("size_t size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type)) my_print("my_type size does not match file", mpi_info.rank);
  if(my_sz !=sizeof(my_type) ||size_sz !=sizeof(size_t)) return 1;
  
  file.read((char*) &verf, sizeof(my_type));
  file.read((char*) &tmp_vers, sizeof(char)*15);

  if(verf != io_verify){
  //equality even though floats as should be identical
    my_print("File read error", mpi_info.rank);
    if(strcmp(tmp_vers, VERSION) !=0 && strcmp(tmp_vers, "IDL data write")!= 0 ) my_print("Incompatible code versions", mpi_info.rank);
    return 1;
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
  if(n_dims_in !=n_dims){
    my_print("Dimensions do not match, aborting read", mpi_info.rank);
    return 1;
  }
  size_t tot_els =1;
  for(size_t i=0;i<n_dims_in;i++){
    file.read((char*) &dim_tmp, sizeof(size_t));
    if(dim_tmp !=dims[i]){
      my_print("Dimensions do not match, aborting read", mpi_info.rank);
      return 1;
    }
    tot_els *= dim_tmp;
    //Check dims
  }
  file.read((char*) &next_block, sizeof(size_t));

  file.read((char *) data , sizeof(my_type)*tot_els);

  return 0;

}

bool my_array::resize(size_t dim, size_t sz){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Note due to using 1-d memory layout both cases require copying all data and therefore briefly memory to store the old and new arrays. However shinking the last dimension does not necessarily require a copy.
*/

  my_print("Attempting to resize", mpi_info.rank);

  if(sz == 0 || sz > MAX_SIZE)return 1;
  //size errors. Don't allow to resize to 0!!!
  if(dim > this->n_dims) return 1;
  if(sz == dims[dim]){
     my_print("Size matches", mpi_info.rank);
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
    for(size_t i=0; i<dim; ++i) part_sz*= dims[i];
    for(size_t i=dim+1; i<n_dims; ++i) part_sz*= dims[i];
    new_data=(my_type*)calloc(part_sz*sz,sizeof(my_type));
    size_t els_to_copy = 1, n_segments = 1;

    // Now we know dim is zeroth or a middle dimension. (1 for n_dims==3, 1 or 2 for n_dims==4.So we copy in chunks
    for(size_t i=0; i<dim; ++i) els_to_copy *= dims[i];
    size_t chunk_sz = els_to_copy;
    (sz> dims[dim])? els_to_copy *= dims[dim] : els_to_copy *= sz;
    for(size_t i=dim+1; i< n_dims; ++i) n_segments *= dims[i];
    for(size_t i=0; i< n_segments; ++i) std::copy(data + i*chunk_sz*dims[dim],data + i*chunk_sz*dims[dim]+ els_to_copy, new_data + i*chunk_sz*sz);

    free(data);
  }

  data = new_data;
  dims[dim] = sz;

  my_print("New size of dim "+mk_str(dim)+  " is " + mk_str(dims[dim]), mpi_info.rank);

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
  size_t part_sz, sub_sz = 1;

/*  if(dim == n_dims -1){
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
/** Find minimum value of data, allows linear search through contiguous memory*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  return *(std::max_element(data+offset, data+total_size));


}
my_type my_array::maxval(std::vector<size_t> &ind, size_t offset){
/** Find minimum value of data, allows linear search through contiguous memory*/

  size_t total_size=get_total_elements();
  if(offset > total_size) return std::numeric_limits<my_type>::max();

  auto it = std::max_element(data+offset, data+total_size);
  ind = get_index_from_offset(it - data);
//  ind.push_back(it - data);

  return *(it);

}


void data_array::construct(){
/** \brief Common parts for all constructors
*
*
*/
  axes=nullptr;
  ax_defined = false;
  time[0]=0; time[1]=1;
  space[0]=0; space[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  
}

void data_array::alloc_ax(const size_t els){
/* \brief Allocate axis memory
*
* Alocate memory for axes and set the defined flag
*/
  if(els > 0 && els <= this->n_dims*MAX_SIZE){
    axes=(my_type*)calloc((els),sizeof(my_type));
    if(axes) ax_defined=true;
  }else{
    my_print("Array size exceeds max. Axes alloc failed", mpi_info.rank);
  }
}

data_array::data_array() : my_array(){
  construct();
}

data_array::data_array(size_t nx, size_t ny, size_t nz, size_t nt) : my_array(nx, ny, nz, nt){
/**Adds axes to a normal rectangular my array*/
  construct();
  size_t els= this->get_total_axis_elements();
  //by now this is setup to work
  alloc_ax(els);
}

data_array::data_array(size_t n_dims, size_t * dims) : my_array(n_dims, dims){

  construct();
  size_t els= this->get_total_axis_elements();
  //by now this is setup to work
  alloc_ax(els);


}

data_array::data_array(std::string filename, bool no_version_check){
/**\brief Create data array from file
*
* Create a data array by reading from the named file. If the file does not exist no memory is allocated. Otherwise it reads the dimensions and sets itself up accordingly \todo Update format
*/

  std::fstream infile;
  infile.open(filename, std::ios::in|std::ios::binary);
  if(!infile.is_open()) return;
  
  char id_in[ID_SIZE];
  char tmp_vers[15];
  my_type verf=0.0;
  int n_dims_in, dim_tmp;

  infile.read(id_in, sizeof(char)*ID_SIZE);
  
  infile.read((char*) &verf, sizeof(my_type));

  infile.read((char*) &tmp_vers, sizeof(char)*15);

  if(verf != io_verify){
  //equality even though floats as should be identical
    my_print("Bugger, file read error", mpi_info.rank);
    if(tmp_vers !=VERSION && strcmp(tmp_vers, "IDL data write")!= 0) my_print("Incompatible code versions", mpi_info.rank);
    return;
  }else{
    if(!no_version_check && tmp_vers !=VERSION) my_print("WARNING: A different code version was used to write this data. Proceed with caution. Fields may not align correctly.", mpi_info.rank);
  }


  infile.read((char*) &n_dims_in, sizeof(int));

  //Now we have the dimensions, construct
  size_t total_data=1, total_axes=0;
  if(n_dims_in >0){
    
    construct();

    this->n_dims = n_dims_in;
    this->dims = (size_t*)malloc(n_dims*sizeof(size_t));
    for(size_t i=0;i<n_dims;i++){
      infile.read((char*) &dim_tmp, sizeof(int));

      dims[i] = (size_t) dim_tmp;
      total_axes += dims[i];
      total_data *= dims[i];
      
    }
    
    data=(my_type*)calloc(total_data,sizeof(my_type));
    if(data) defined=true;
    axes=(my_type*)calloc(total_axes,sizeof(my_type));
    if(axes) ax_defined=true;

    //Finally read in data and axes using normal routines
    infile.seekg(0, std::ios::beg);
    this->read_from_file(infile, no_version_check);

  
  }else{
    my_print("Invalid dimensionality in input file", mpi_info.rank);
  
  }


}

data_array::~data_array(){
/**Similarly destructor automatically calls destructor for my_array and frees axes*/
  if(axes) free(axes);
  axes = nullptr; // technically unnecessary as desctructor deletes members.
}

data_array & data_array::operator=(const data_array& src){
  
  my_array::operator=(std::move(src));
  if(this->axes) free(axes);
  this->construct();

  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
  //Allocate axis memory
  if(axes) std::copy(src.axes, src.axes+els, this->axes);
  //Copy axes
  copy_ids(src);
  return *this;

}

data_array::data_array(const data_array &src) : my_array(src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data \todo move constructor
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

my_type * data_array::get_axis(size_t dim, size_t & length){
/**  \brief Get pointer to axis
*
*Returns pointer to given axis and its length. If axes don't exist or dimension is out of range, returns nullptr
*/

  if(!ax_defined || (dim >= n_dims)) return nullptr;

  long index = get_axis_index(dim, 0);
  //Get index of 0th element
  length = get_length(dim);
  if(index != -1) return axes + index;
  else return nullptr;

}

long data_array::get_axis_index(size_t dim, size_t pt){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. Let this function do all bounds checks.
*/

  if(dim >=n_dims || pt >=get_length(dim)) return -1;
  //Out of range error
  
  long offset = 0;
  for(size_t i=0; i< dim; i++) offset +=dims[i];
  return offset + pt;

}

my_type data_array::get_axis_element(size_t dim, size_t pt){
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
  if(axis && len >1) return std::abs(axis[0]-axis[dims[i]-1])/(dims[i]-1);
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

void data_array::make_linear_axis(size_t dim, float res, size_t offset){
/**\brief Make an axis
*
*Generates a linear axis for dimension dim, with resolution res, starting at value of  - offset*res @param dim Dimension to build axis for @param res Axis resolution @param offset Number of grid cells to shift downwards (leftwards) by
*/

  size_t len;
  my_type * ax_ptr = get_axis(dim, len);

  for(size_t i=0; i<len; i++) *(ax_ptr +i) = ((float) (i-offset)) * res;

}

bool data_array::write_to_file(std::fstream &file){
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
  if(file.tellg() != next_location) write_err=1;

  size_t hdr_start = next_location + sizeof(size_t);

  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write(block_id, sizeof(char)*ID_SIZE);

  if(file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & hdr_start, sizeof(size_t));
  //Finish with position of start of footer!
  return 0;

}

bool data_array::write_section_to_file(std::fstream &file, std::vector<my_type> limits){
/** \brief Print section of array to file
*
*Prints the section defined by the vector limits to supplied file. Limits should contain AXIS values. To use one dimension entire supply values less/greater than min and max axis values. See data_array::write_to_file for format
*/

  if(!file.is_open()) return 1;
  bool write_err = 0;

  if(limits.size() != 2*n_dims){
    my_print("Limits vector size does not match array!", mpi_info.rank);
    return 1;
  }

  //Identify limits of segment from axes
  std::vector<size_t> index_limits = this->get_bounds(limits);

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

  if(file.tellg() != next_location) write_err=1;

  size_t hdr_start = next_location  + sizeof(size_t);
  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  std::cout<<block_id<<'\n';
  file.write(block_id, sizeof(char)*ID_SIZE);

  if(file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & hdr_start, sizeof(size_t));

  return 0;
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
  size_t where_val;
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
*Read a file into a pre-sized data array See also data_array::data_array(std::string filename, bool no_version_check);
*/
/*
*First write the my_array section, see my_array::write_to_file then add the following
* Next_block axes
* Next_block Block_id Prev_block
*/

  bool err;

  if(file.good()) err=my_array::read_from_file(file, no_version_check);
  else my_print("Read error!", mpi_info.rank);
  if(err){
    my_print("File read failed", mpi_info.rank);
    return err;
  }
  //call parent class to read data, checking we read id ok first

  size_t next_block=0, end_block=0;
  file.read((char*) &next_block, sizeof(size_t));

  size_t tot_els = get_total_axis_elements();

  if(!err) file.read((char *) this->axes , sizeof(my_type)*(tot_els));
  //If we managed to get dims etc, continue on to get axes

  file.read((char*) &next_block, sizeof(size_t));

  file.seekg(-1*sizeof(size_t), file.end);
  file.read((char*) &end_block, sizeof(size_t));
  //First read the block ID
  char id_in[ID_SIZE];
  file.seekg(end_block);
  if(file) file.read(id_in, sizeof(char)*ID_SIZE);
  strcpy(this->block_id, id_in);

  return err;

}

bool data_array::fft_me(data_array * data_out){
/** \brief FFT data_array
*
* Data and axes in this object are FFT'd using FFTW and stored into the instance pointed to by data_out. Data_out must be created with correct dimensions first, but we check and return error (1) if it is not so.
*/

  if(!data_out->is_good()){
    my_print("Output array for FFT undefined", mpi_info.rank);
    return 1;
  }
  if(data_out->n_dims != this->n_dims){
    my_print("Wrong output dimensions for FFT", mpi_info.rank);
    return 1;
  }
  for(size_t i=0; i<n_dims;++i){
    if(data_out->dims[i] != this->dims[i]){
      my_print("Wrong output dimensions for FFT", mpi_info.rank);
      return 1;
    }
  }

  data_out->copy_ids(*this);
  size_t total_size=1; /* Total number of elements in array*/
  for(size_t i=0; i<n_dims;++i) total_size *= dims[i];

  int fft_dim ;
  fft_dim = 1;/* Dimension to FFT over, if required*/

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
  
  bool err=false;
  err = data_out->populate_mirror_fastest(result, total_size);
  //Copy result into out array

  if(err){
    my_print("Error populating result array", mpi_info.rank);
    return 1;
  }

  ADD_FFTW(destroy_plan)(p);
  ADD_FFTW(free)(in);
  ADD_FFTW(free)(out);
  ADD_FFTW(free)(result);
  //Destroy stuff we don't need

  //do shifts for all but 0th dim
  size_t shft = 0;
  for(size_t i=1; i< n_dims; i++){
    shft = data_out->get_dims(i)/2;
    data_out->shift(i, shft, 0);
  }


  my_type * tmp_axis;
  float N2, res;
  size_t len;

  for(size_t i=0;i<n_dims;++i){
  //Loop over the dimensions and construct each axis in place. We KNOW now that data_out has correct dimensions. We checked before getting here. We don't need to check len. It is the same as dims[i] each time. Perhaps we will anyway? :)
    tmp_axis = data_out->get_axis(i, len);
    if(len != dims[i]) return 1;

    N2 = ((float) dims[i])/2.0;
    res = this->get_res(i);
    for(size_t j= 0; j< dims[i]; j++) *(tmp_axis + j) = pi * ((float)j -N2)/N2/res;
  }

  return 0;

}

bool data_array::populate_mirror_fastest(my_type * result_in, size_t total_els){

//  int total_size=1;
//  for(int i=0; i< n_dims; i++) total_size*=dims[i];

//  return this->populate_data(result_in, total_els);
  size_t last_size = dims[0]/2 + 1;
  size_t num_strides = total_els/dims[0];
  //std::cout<<total_els<<" "<<last_size<<" "<<dims[0]<<" "<<num_strides<<" "<<this->get_total_elements()<<'\n';
 
  for(size_t i=0; i< num_strides; i++){
    //std::cout<<i*dims[0]+last_size-1<<" "<<i*dims[0]+last_size-1+last_size<<'\n';
    std::copy(result_in+ i*last_size, result_in+(i+1)*last_size -1, data+i*dims[0]+last_size-1);
    std::reverse_copy(result_in+ i*last_size, result_in+(i+1)*last_size-1, data+i*dims[0]+1);
  }
  return 0;
}

void data_array::copy_ids( const data_array &src){
/** Copies ID fields from src array to this */

  strcpy(this->block_id, src.block_id);
  for(size_t i=0; i < 2; ++i) this->time[i] = src.time[i];
  //  std::copy(src->time, src->time + 2, this->time);
  for(size_t i=0; i < 2; ++i) this->space[i] = src.space[i];
}

bool data_array::check_ids( const data_array & src){
/** Checks ID fields match src */

  bool err=false;
  if(strcmp(this->block_id, src.block_id) != 0) err =true;
  for(size_t i=0; i< 2; i++) if(src.time[i] != this->time[i]) err=true;
  for(size_t i=0; i < 2; ++i) if(this->space[i] != src.space[i]) err=true;

  return err;
}

bool data_array::resize(size_t dim, size_t sz){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Similarly for axis elements. See my_array::resize() for more.
*/

  size_t old_sz = this->get_dims(dim);
  //call my_array::resize to resize data...
  bool err = my_array::resize(dim, sz);

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

