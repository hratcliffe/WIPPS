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
#include <fftw3.h>
#include <cmath>
#include "support.h"
#include "my_array.h"

extern const mpi_info_struc mpi_info;


/*void my_array::tmp_function(){

  std::cout<< mpi_info.rank<<std::endl;
    mpi_info.rank +=1;
    std::cout<< mpi_info.rank<<std::endl;
  if(mpi_info.rank ==0) std::cout<< &mpi_info<<std::endl;

}*/


my_array::my_array(){
/** Default constructor*/
  construct();
}

void my_array::construct(){
/** \brief Shared contructor code
*
*Sets default values of things in case of constructor omissions
*/
  defined = false;
  ragged=false;
  n_dims = 0;
  data=nullptr;

}

my_array::my_array(int nx, int ny, int nz, int nt){
/** \brief 2-d rectangular array
*
*Sets up a n-d rectangular array of nx*ny and allocates data arrays. If either of first 2 sizes is zero or exceeds MAX_SIZE, print error and exit. Otherwise construct 2, 3, 4d as dims specify
*/
  construct();
  
  if(ny==0 || nx==0){
    my_print("Array cannot have 0 dim", mpi_info.rank);
    return;
  }
  if(ny> MAX_SIZE || nx>MAX_SIZE|| nz>MAX_SIZE || nt>MAX_SIZE){
    my_print("Array size exceeds MAX_SIZE of "+mk_str(MAX_SIZE), mpi_info.rank);
    return;
  }
  if(nz==0 && nt==0){
    n_dims = 2;
    this->dims = (int*)malloc(n_dims*sizeof(int));
    dims[0]=nx;
    dims[1]=ny;

    data=(my_type*)calloc(nx*ny,sizeof(my_type));
  }else if (nt==0){
    n_dims = 3;
    this->dims = (int*)malloc(n_dims*sizeof(int));
    dims[0]=nx;
    dims[1]=ny;
    dims[2]=nz;

    data=(my_type*)calloc(nx*ny*nz,sizeof(my_type));
  
  }else if (nz>0 && nt>0){
    n_dims = 4;
    this->dims = (int*)malloc(n_dims*sizeof(int));
    dims[0]=nx;
    dims[1]=ny;
    dims[2]=nz;
    dims[3]=nt;

    data=(my_type*)calloc(nx*ny*nz*nt,sizeof(my_type));
  
  }else{
    my_print("Array cannot have 0 dim", mpi_info.rank);
    return;

  }
  if(data){
    defined = true;
  }
}

my_array::my_array(int * row_len, int ny){
/** \brief 2-d ragged array
*
*Sets up a 2-d array containing ny rows of lengths given in row_len and allocates data arrays. If ny is zero or exceeds MAX_SIZE, or any row_len is, print error and exit
*/

  construct();
  
  int nx_min=MAX_SIZE, nx_max=0;
  for(int i=0; i< ny; ++i){
    if(row_len[i]< nx_min) nx_min = row_len[i];
    if(row_len[i]> nx_max) nx_max = row_len[i];
  }
  if(ny==0 || nx_min==0){
    my_print("Array size cannot be 0", mpi_info.rank);
    return;
  }
  if(ny> MAX_SIZE || nx_max >MAX_SIZE){
    my_print("Array size exceeds MAX_SIZE of "+mk_str(MAX_SIZE), mpi_info.rank);
    return;
  }

  n_dims = 2;
  dims = (int*)malloc(n_dims*sizeof(int));
  dims[0]=0;
  dims[1]=ny;
  ragged = true;

  this->row_lengths = (int*) malloc(ny*sizeof(int));
  this->cumulative_row_lengths = (int*) malloc(ny*sizeof(int));

  memcpy((void *) this->row_lengths, (void *)row_len, ny*sizeof(int));
  
  cumulative_row_lengths[0] = 0;

  for(int i=1; i<ny;++i) cumulative_row_lengths[i] = cumulative_row_lengths[i-1] + this->row_lengths[i-1];
  
  data=(my_type*)calloc((cumulative_row_lengths[ny-1] + row_lengths[ny-1]), sizeof(my_type));

  if(data){
    defined = true;
  }
}

my_array::~my_array(){
/** Clean up explicit allocations
*
*
*/
  if(data) free(data);
  data = NULL; // technically unnecessary as desctructor deletes memebers. 
  if(dims) free(dims);
  if(ragged){
    if(cumulative_row_lengths) free(cumulative_row_lengths);
    if(row_lengths) free(row_lengths);
  }
}

int my_array::get_index(int nx){
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
    my_print("Wrong array dimension, check get_index calls", mpi_info.rank);
    return -1;

  }
  if((nx < dims[0])){
    return nx;
  }else{
    return -1;
  }

}
int my_array::get_index(int nx, int ny){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args

* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>
Etc
*/

  if(n_dims != 2){
    my_print("Wrong array dimension, check get_index calls", mpi_info.rank);
    return -1;

  }
  if(!ragged){
    if((nx < dims[0]) && (ny<dims[1])){
      return ny*dims[0] + nx;
    }else{
      return -1;
    }
  }else{
    //have to check specific row length...
    if((ny<dims[1]) && (nx < row_lengths[ny])){
      return cumulative_row_lengths[ny] + nx;
    }else{
      return -1;
    }
  
  }
}
int my_array::get_index(int nx, int ny, int nz){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args
*/

  if(n_dims != 3){
    return -1;
    my_print("Wrong array dimension, check get_index calls", mpi_info.rank);
  }

  if(!ragged){
    if((nx < dims[0]) && (ny<dims[1]) && ((nz<dims[2]))){
      return (nz*dims[1]+ ny)*dims[0] + nx;
    }else{
      return -1;
    }
  }else{
    return -1;
     /**\todo Can we do this case elegantly??*/
  }
}
int my_array::get_index(int nx, int ny, int nz, int nt){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. NB. Let this function do all bounds checks. Just call it plain. This function is called often so we make it as simple as possible and write one for each number of args
*/

  if(n_dims != 4){
    return -1;
    my_print("Wrong array dimension, check get_index calls", mpi_info.rank);
  }
  if(!ragged){
    if((nx < dims[0]) && (ny<dims[1])&& ((nz<dims[2]))&& ((nt<dims[3]))){
      return ((nt*dims[2] +nz)*dims[1]+ ny)*dims[0] + nx;
    }else{
      return -1;
    }
  }else{
     /**\todo Can we do this case elegantly??*/
    return -1;

  }
}

int my_array::get_dims(){
  return n_dims;
}
int my_array::get_dims(int dim){
/** \brief Return dims[dim]
*
*WARNING Do NOT use this to get sizes of a ragged array. It won't work. Use get_length instead
*/
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

int my_array::get_length(int dim){
/** \brief Get size of dimension dim
*
*For a rectangular array returns the stored dim, else returns the row length for ragged
*/
  if(!ragged) return get_dims(dim);
  else{
    if(n_dims==2 && dim < dims[1]) return row_lengths[dim];
    else return 0;
  }

}

my_type my_array::get_element(int nx){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  int ind = get_index(nx);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(int nx, int ny){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  int ind = get_index(nx, ny);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }
}
my_type my_array::get_element(int nx, int ny, int nz){
/** Return element at nx, ny, nz. Out of range etc will return 0.0*/
  int ind = get_index(nx, ny, nz);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}
my_type my_array::get_element(int nx, int ny, int nz, int nt){
/** Return element at nx, ny, nz, nt. Out of range etc will return 0.0*/
  int ind = get_index(nx, ny, nz, nt);
  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }

}

int my_array::get_total_elements(){
/** Return total size of array */
  int tot_els=1;

  if(!ragged){
    for(int i=0; i<n_dims;++i) tot_els *=dims[i];
  }else{
    tot_els = cumulative_row_lengths[dims[n_dims-1]-1]+row_lengths[dims[n_dims-1]-1];
  }

  return tot_els;

}

bool my_array::set_element(int nx, int ny, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, @return 1 if out of range, wrong number of args, 0 else.
*/

  int index = get_index(nx, ny);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}
bool my_array::set_element(int nx, int ny, int nz, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, @return 1 if out of range, wrong number of args, 0 else.
*/

  int index = get_index(nx, ny, nz);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}
bool my_array::set_element(int nx, int ny, int nz, int nt, my_type val){
/** \brief Sets array element
*
*Sets elements at nx, ny, nz, nt, @return 1 if out of range, wrong number of args, 0 else.
*/

  int index = get_index(nx, ny, nz, nt);
  if(index >= 0){
    data[index] = val;
    return 0;
  }else{
    return 1;
  }
  
}

bool my_array::populate_data(my_type * dat_in, int n_tot){
/** \brief Fill array
*
*Populates data with the total number of elements specified, as long as n_tot is less than the product of dims. Parameter needed so we can't overflow dat_int. Assumes same row-column order etc etc @return 0 (sucess) 1 (error)
*/

  int tot_els = get_total_elements();
//  if(n_tot > tot_els) return 1;
  if(n_tot < tot_els) tot_els = n_tot;
  //Use min of n_tot, tot_els

  void * tmp = (void *) this->data;
  if(!tmp) return 1;

  memcpy (tmp , dat_in, tot_els*sizeof(my_type));

  return 0;

}

bool my_array::populate_row(void * dat_in, int nx, int y_row){
/** /brief Fill row
*
* Fills the row y_row of array from dat_in to length of nx. Nx is param for sanity so we can't overflow dat_in. @return 0 (sucess) 1 (error) \todo Extend to larger dimensions
*/

  if(nx != dims[0]) return 1;
  if(y_row > dims[1] || y_row < 0) return 1;
  int indx = get_index(0, y_row);

  void * tmp = (void *) (data+indx);

  if(indx==-1 || !tmp) return 1;
  memcpy (tmp , dat_in, nx*sizeof(my_type));

  return 0;
}

bool my_array::write_to_file(std::fstream &file){
/**Takes the version etc info then the whole data array and writes as a stream. It's not portable to other machines necessarily due to float sizes and endianness. It'll do. We'll start the file with a known float for confirmation.
  *
  *IMPORTANT: this VERSION specifier links output files to code. If modifying output or order commit and clean build before using. @return 0 (sucess) 1 (error)
*/
  if(!file.is_open() || (this->data ==nullptr)) return 1;
  const char tmp_vers[15] = VERSION;

  file.write((char*) &io_verify, sizeof(my_type));
  file.write((char*) &tmp_vers, sizeof(char)*15);

  //Code version...

  int total_size = get_total_elements();
  //dimension info
  if(!ragged){
    file.write((char*) &n_dims, sizeof(int));
    int dim_tmp;
    for(int i=0;i<n_dims;i++){
      dim_tmp = dims[i];
      file.write((char*) &dim_tmp, sizeof(int));
    }
  }else{
    //do different if we have ragged array...
    int n_dims_new = -1*n_dims;
    file.write((char*) &n_dims_new, sizeof(int));
    int dim_tmp;
      for(int i=0;i<n_dims;i++){
      dim_tmp = dims[i];
      file.write((char*) &dim_tmp, sizeof(int));
    }

    for(int i=0;i<dims[n_dims-1];i++){
      dim_tmp = row_lengths[i];
      file.write((char*) &dim_tmp, sizeof(int));
    }

  }
  file.write((char *) data , sizeof(my_type)*total_size);

  return 0;

}


bool my_array::write_section_to_file(std::fstream &file, std::vector<int> bounds){
/**Takes the version etc info then the section of the data array and writes as a stream. It's not portable to other machines necessarily due to float sizes and endianness. It'll do. We'll start the file with a known float for confirmation.
  * We use lazy method of get_element for each element, less prone to offset errors and memory is already much faster than disk
  *
  *IMPORTANT: this VERSION specifier links output files to code. If modifying output or order commit and clean build before using. @return 0 (sucess) 1 (error)
*/

  if(!file.is_open() || (this->data ==nullptr)) return 1;
  if(bounds.size() != n_dims*2) return 1;
  for(int i=0; i< n_dims; i++){
    if(bounds[2*i] < 0 || bounds[2*i+1] > dims[i]) return 1;
  }
  
  const char tmp_vers[15] = VERSION;

  file.write((char*) &io_verify, sizeof(my_type));
  file.write((char*) &tmp_vers, sizeof(char)*15);
  //Code version...

  int total_size = get_total_elements();
  //dimension info
  if(!ragged){
    file.write((char*) &n_dims, sizeof(int));
    int dim_tmp;
    for(int i=0;i<n_dims;i++){
      dim_tmp = bounds[i*2+1]-bounds[i*2];
      file.write((char*) &dim_tmp, sizeof(int));
    }
  }else{
    //do different if we have ragged array...
    int n_dims_new = -1*n_dims;
    file.write((char*) &n_dims_new, sizeof(int));
    int dim_tmp;
      for(int i=0;i<n_dims;i++){
      dim_tmp = dims[i];
      file.write((char*) &dim_tmp, sizeof(int));
    }

    for(int i=0;i<dims[n_dims-1];i++){
      dim_tmp = row_lengths[i];
      file.write((char*) &dim_tmp, sizeof(int));
    }

  }
  my_type element;
  if(n_dims ==1){
    for(int i= bounds[0]; i< bounds[1]; i++){
      element = get_element(i);
      file.write((char *) &element, sizeof(my_type));
    }
  }else if(n_dims ==2){
    for(int j= bounds[2]; j< bounds[3]; j++){
      for(int i= bounds[0]; i< bounds[1]; i++){
        element = get_element(i, j);
        file.write((char *)  &element, sizeof(my_type));
      }
    }
  }else if(n_dims ==3){
    for(int k= bounds[4]; k< bounds[5]; k++){
      for(int j= bounds[2]; j< bounds[3]; j++){
        for(int i= bounds[0]; i< bounds[1]; i++){
          element  = get_element(i, j, k);
          file.write((char *) &element, sizeof(my_type));
        }
      }
    }
  }else if(n_dims ==4){
    for(int l= bounds[6]; l< bounds[7]; l++){
      for(int k= bounds[4]; k< bounds[5]; k++){
        for(int j= bounds[2]; j< bounds[3]; j++){
          for(int i= bounds[0]; i< bounds[1]; i++){
            element = get_element(i, j, k, l);
            file.write((char *) &element , sizeof(my_type));
          }
        }
      }
    }
  }
  return 0;

}


bool my_array::read_from_file(std::fstream &file, bool no_version_check){
/** \brief File read
*
*Reads block id, data and axes from a file. Requires the dimensions of array to be already setup and will check for consistency with those in file
*/

  char tmp_vers[15];
  my_type verf=0.0;

  int n_dims_in, dim_tmp;

  file.read((char*) &verf, sizeof(my_type));
  file.read((char*) &tmp_vers, sizeof(char)*15);

  if(verf != io_verify){
  //equality even though floats as should be identical
    my_print("Bugger, file read error", mpi_info.rank);
    if(tmp_vers !=VERSION && strcmp(tmp_vers, "IDL data write")!= 0 ) my_print("Incompatible code versions", mpi_info.rank);
    return 1;
  }else{
    if(!no_version_check && tmp_vers !=VERSION) my_print("WARNING: A different code version was used to write this data. Proceed with caution. Fields may not align correctly.", mpi_info.rank);
  }

  file.read((char*) &n_dims_in, sizeof(int));
  if(n_dims_in !=n_dims){
    my_print("Dimensions do not match, aborting read", mpi_info.rank);
    return 1;
  }
  int tot_els =1;
  for(int i=0;i<n_dims_in;i++){
    file.read((char*) &dim_tmp, sizeof(int));
    if(dim_tmp !=dims[i]){
      my_print("Dimensions do not match, aborting read", mpi_info.rank);
      return 1;
    }
    tot_els *= dims[i];
  }
  
  file.read((char *) data , sizeof(my_type)*tot_els);

  return 0;

}

bool my_array::resize(int dim, int sz){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Note due to using 1-d memory layout both cases require copying all data and therefore briefly memory to store the old and new arrays. However shinking the last dimension does not necessarily require a copy. Note cannot be called on ragged array. NOTE dim runs from 1 to number of dims
*/

  my_print("Attempting to resize", mpi_info.rank);

  if(sz < 0 || sz > MAX_SIZE) return 1;
  //size errors
  if(dim > this->n_dims) return 1;
  if(dim>=0 && sz == dims[dim]){
     my_print("Size matches", mpi_info.rank);
     return 1;
  }
  if(ragged){
   my_print("Cannot resize a ragged array. Create new and copy", mpi_info.rank);
     return 1;
  }
  
  my_type * new_data;
  int part_sz = 1;

  if(dim == n_dims-1){
    //special case as we can shrink and maybe grow without copy
    for(int i=0; i<n_dims-1; ++i) part_sz*= dims[i];
    //product of all other dims

    new_data = (my_type *) realloc((void*) this->data, part_sz*sz*sizeof(my_type));
    if(!new_data){
      my_print("Failed to reallocate memory", mpi_info.rank);
      return 1;
      //failure. leave as was.
    }
    int new_els = part_sz*(sz - dims[n_dims-1]);
    
    if(new_els > 0) memset((void*)(new_data + part_sz*dims[n_dims-1]), 0.0, new_els*sizeof(my_type));
    //zero new elements
  }
  else{
    //have to allocate a new block and copy across.
    for(int i=0; i<dim-1; ++i) part_sz*= dims[i];
    for(int i=dim+1; i<n_dims; ++i) part_sz*= dims[i];

    new_data=(my_type*)calloc(part_sz*sz,sizeof(my_type));
    int els_to_copy = 1, n_segments = 1;

    // Now we know dim is zeroth or a middle dimension. (1 for n_dims==3, 1 or 2 for n_dims==4.So we copy in chunks
    for(int i=0; i<dim; ++i) els_to_copy *= dims[i];
    int chunk_sz = els_to_copy;
    (sz> dims[dim])? els_to_copy *= dims[dim] : els_to_copy *= sz;
    for(int i=dim+1; i< n_dims; ++i) n_segments *= dims[i];
    for(int i=0; i< n_segments; ++i) memcpy((void*)(data + i*chunk_sz*dims[dim]), (void*)(new_data + i*chunk_sz*sz), els_to_copy);

    free(data);
  }

  data = new_data;
  dims[dim] = sz;

  my_print("New size of dim "+mk_str(dim)+  " is " + mk_str(dims[dim]), mpi_info.rank);

  return 0;

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

data_array::data_array(int nx, int ny) : my_array(nx,ny){
/**Adds axes to a normal rectangular my array*/

  construct();
  axes=(my_type*)calloc((nx+ny),sizeof(my_type));
  if(axes) ax_defined=true;

}
data_array::data_array(int nx, int ny, int nz) : my_array(nx,ny, nz){
/**Adds axes to a normal rectangular my array*/

  construct();
  axes=(my_type*)calloc((nx+ny+nz),sizeof(my_type));
  if(axes) ax_defined=true;

}
data_array::data_array(int nx, int ny, int nz, int nt) : my_array(nx,ny, nz, nt){
/**Adds axes to a normal rectangular my array*/

  construct();
  axes=(my_type*)calloc((nx+ny+nz+nt),sizeof(my_type));
  if(axes) ax_defined=true;

}

data_array::data_array(int * row_lengths, int ny): my_array(row_lengths,ny){
/** Adds axes to a ragged my_array One per row in this case...*/

  construct();

  int tot_els = cumulative_row_lengths[dims[n_dims-1]-1]+row_lengths[dims[n_dims-1]-1];

  axes=(my_type*)calloc(tot_els, sizeof(my_type));
  if(axes) ax_defined=true;

}

data_array::data_array(std::string filename, bool no_version_check){
/**\brief Create data array from file
*
* Create a data array by reading from the named file. If the file does not exist no memory is allocated. Otherwise it reads the dimensions and sets itself up accordingly
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
  int total_data=1, total_axes=0;
  if(n_dims_in >0){
    
    construct();

    this->n_dims = n_dims_in;
    this->dims = (int*)malloc(n_dims*sizeof(int));
    for(int i=0;i<n_dims;i++){
      infile.read((char*) &dim_tmp, sizeof(int));

      dims[i] = dim_tmp;
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

  
  }else if(n_dims_in < 0){
    //This means ragged array
    my_print("I didn't finish this because i am useless!!!!!!!!!!!!!!!!!!!!!!!!!!!", mpi_info.rank);
  
  }else{
    my_print("Invalid dimensionality in input file", mpi_info.rank);
  
  }


}

data_array::~data_array(){
/**Similarly destructor automatically calls destructor for my_array and frees axes*/

  if(axes) free(axes);
  axes = NULL; // technically unnecessary as desctructor deletes members.
  

}

my_type * data_array::get_axis(int dim, int & length){
/**  \brief Get pointer to axis
*
*Returns pointer to given axis and its length. If axes don't exist or dimension is out of range, returns nullptr
*/

  if(!ax_defined || (dim >= n_dims)) return nullptr;

  int index = get_axis_index(dim, 0);
  //Get index of 0th element
  length = get_length(dim);
  if(index != -1) return axes + index;
  else return nullptr;

}

int data_array::get_axis_index(int dim, int pt){
/** \brief Get index for location
*
*Takes care of all bounds checking and disposition in memory. Returns -1 if out of range of any sort, otherwise, suitable index. Let this function do all bounds checks.
*/

  if(dim < 0 || dim >=n_dims || pt >=get_length(dim)) return -1;
  //Out of range error
  
  int offset = 0;
  if(!ragged){
  // Rectangular, skip over other dims
    for(int i=0; i< dim; i++) offset +=dims[i];
  
  }else{
  // Ragged, skip axes for other rows
    offset = cumulative_row_lengths[dim];
  }
  
  return offset + pt;

}

my_type data_array::get_axis_element(int dim, int pt){
/** \brief Get axis value
*
*Returns value at pt in dimension dim if in range, else 0.0
*/

  int ind = get_axis_index(dim, pt);
  if(ind  != -1){
    return axes[ind];
  }else{
    return 0.0;
  }

}

bool data_array::set_axis_element(int dim, int pt, my_type val){
/** \brief Sets array element
*
*Sets elements at pt on dimension dim, @return 1 if out of range, 0 else.
*/

  int index = get_axis_index(dim, pt);
  if(index >= 0){
    axes[index] = val;
    return 0;
  }else{
    return 1;
  }

}

bool data_array::populate_axis(int dim, my_type * dat_in, int n_tot){
/** \brief Fill axis
*
*Populates axis with the total number of elements specified, as long as n_tot is less than the specified dim. Parameter needed so we can't overflow dat_in. @return 0 (sucess) 1 (error)
*/

  if(dim < 0 || dim >=n_dims) return -1;
  //Out of range error

  int tot_els = dims[dim];
  if(n_tot < tot_els) tot_els = n_tot;
  //Min of n_tot and tot_els

  void * tmp = (void *) this->axes;
  if(!tmp) return 1;

  memcpy (tmp , dat_in, tot_els*sizeof(my_type));

  return 0;

}   
float data_array::get_res(int i){
/**Return resolution of axis on dimension i. Assumes linear etc etc. If axis is undefined or zero or one in length, return 1.0 */
  int len;
  my_type * axis = this->get_axis(i, len);
  if(axis && len >1) return std::abs(axis[0]-axis[dims[i]-1])/(dims[i]-1);
  else return 1.0;

}

int data_array::get_total_axis_elements(){
/** \brief Return total axes length
*
*Sums number of total elements in all axes
*/
  int tot_els=0;

  if(!ragged){
    for(int i=0; i<n_dims;++i) tot_els +=dims[i];
  }else{
    tot_els = cumulative_row_lengths[dims[n_dims-1]-1]+row_lengths[dims[n_dims-1]-1];
  }

  return tot_els;

}

void data_array::make_linear_axis(int dim, float res, int offset){
/**\brief Make an axis
*
*Generates a linear axis for dimension dim, with resolution res, starting at value of  - offset*res @param dim Dimension to build axis for @param res Axis resolution @param offset Number of grid cells to shift downwards (leftwards) by
*/

  int len;
  my_type * ax_ptr = get_axis(dim, len);

  for(int i=0; i<len; i++) *(ax_ptr +i) = ((float) (i-offset)) * res;

}

bool data_array::write_to_file(std::fstream &file){

/**IMPORTANT: the VERSION specifier links output files to code. If modifying output or order commit and clean build before using.
*/

  if(!file.is_open()) return 1;
  file.write(block_id, sizeof(char)*ID_SIZE);

  my_array::write_to_file(file);
  //call base class method to write that data.

  file.write((char *) axes , sizeof(my_type)*(get_total_axis_elements()));
  //Add axes.

  return 0;

}

bool data_array::write_section_to_file(std::fstream &file, std::vector<my_type> limits){
/** \brief Print section of array to file
*
*Prints the section defined by the vector limits to supplied file. Limits should contain AXIS values. To use one dimension entire supply values less/greater than min and max axis values.
*/

  if(!file.is_open()) return 1;
  if(limits.size() != 2*n_dims){
    my_print("Limits vector size does not match array!", mpi_info.rank);
    return 1;
  }
  file.write(block_id, sizeof(char)*ID_SIZE);

  //Identify limits of segment from axes
  std::vector<int> index_limits = this->get_bounds(limits);


  my_array::write_section_to_file(file, index_limits);
  //call base class method to write that data.

//  file.write((char *) axes ,sizeof(my_type)*(get_total_axis_elements()));
  int len;
  for(int i=0; i< n_dims; i++){
    file.write((char *) get_axis(i, len)+index_limits[2*i], sizeof(my_type)*(index_limits[2*i +1]-index_limits[2*i]));

  }
  int i=1;
  std::cout<<index_limits[2*i]<<" "<<index_limits[2*i +1]<<'\n';
  for(int j=index_limits[2*i]; j<index_limits[2*i +1]; j++ ) std::cout<<*(get_axis(i, len)+j)<<" ";

  //Add axes.

  return 0;
}

std::vector<int> data_array::get_bounds(std::vector<my_type> limits){

  std::vector<int> index_limits;

  if(limits.size() != 2*n_dims){
    my_print("Limits vector size does not match array!", mpi_info.rank);
    return index_limits;
  }

  index_limits.resize(n_dims*2);

  int len, index;
  my_type * ax_start;
  int where_val;
  for(int i=0; i< n_dims; i++){
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
/** \brief Test file read
*
*Spams stuff to screen to check read/write
*/

  bool err;
  //First read the block ID
  char id_in[ID_SIZE];

  file.read(id_in, sizeof(char)*ID_SIZE);
  strcpy(this->block_id, id_in);

  if(file.good())err=my_array::read_from_file(file, no_version_check);
  else my_print("Read error!", mpi_info.rank);
  //call parent class to read data, checking we read id ok first

  int tot_els = 0;
  for(int i=0;i<n_dims;i++) tot_els += dims[i];

  if(!err) file.read((char *) this->axes , sizeof(my_type)*(tot_els));
  //If we managed to get dims etc, continue on to get axes

  return err;

}

bool data_array::fft_me(data_array * data_out){
/** \brief FFT data_array
*
* Data and axes in this object are FFT'd using FFTW and stored into the instance pointed to by data_out. Data_out must be created with correct dimensions first, but we check and return error (1) if it is not so. \todo Add 3, 4 dimensions \todo Check handling of odd vs even total sizes
*/

  if(!data_out->is_good()){
    my_print("Output array for FFT undefined", mpi_info.rank);
    return 1;
  }
  if(data_out->n_dims != this->n_dims){
    my_print("Wrong output dimensions for FFT", mpi_info.rank);
    return 1;
  }
  for(int i=0; i<n_dims;++i){
    if(data_out->dims[i] != this->dims[i]){
      my_print("Wrong output dimensions for FFT", mpi_info.rank);
      return 1;
    }
  }

  data_out->copy_ids(this);
  int total_size=1; /* Total number of elements in array*/
  for(int i=0; i<n_dims;++i) total_size *= dims[i];

  int fft_dim ;
  fft_dim = 1;/* Dimension to FFT over, if required*/

  ADD_FFTW(plan) p;
  cplx_type *out;
  my_type * in, *result;

  in = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);
  //my_type should match the used FFTW library, so no type conversion necessary
  out = (cplx_type *) ADD_FFTW(malloc)(sizeof(cplx_type) * total_size);

  result = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);

  /** \todo Possibly this bit can be genericised?*/

  p = ADD_FFTW(plan_dft_r2c)(n_dims, dims, in, out, FFTW_ESTIMATE);

/*  if(n_dims == 1 || (n_dims == 2 && dims[1] == 1) ){
    p = ADD_FFTW(plan_dft_r2c_1d)(dims[0], in, out, FFTW_ESTIMATE);

  }else if(n_dims == 2){
    p = ADD_FFTW(plan_dft_r2c_2d)(dims[0], dims[1], in, out, FFTW_ESTIMATE);

  }else{
    my_print("FFT of more than 2-d arrays not added yet", mpi_info.rank);

    return 1;
  }
*/
  //copy data into in. Because the plan creation changes in, so we don't want to feed our actual data array in, and it's safer to run the plan with the memory block it was created with
  std::copy(this->data, this->data+total_size, in);

  ADD_FFTW(execute)(p);
  //Execute the plan

  cplx_type * addr;
  addr = out;
  //because double indirection is messy and cplx type is currently a 2-element array of floats
  int middle = total_size/2;
  if(middle*2 !=total_size) middle++;
  //odd or even total length
  *(result) = 0.0; //FAKENUMBERS

  for(int i=0; i< middle ; i++){
    *(result+i+middle) = (my_type)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1]));
    *(result-i+middle) = *(result+i+middle);
    //(my_type)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1]));
    
    addr++;
  }
  //Absolute square of out array to produce final result of type my_type
  /** Make sure for 2-d we are making shape we expect!*/
  
  bool err=false;
  err = data_out->populate_data(result, total_size);
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

  my_type * tmp_axis;
  float N2, res;
  int len;

  for(int i=0;i<n_dims;++i){
  //Loop over the dimensions and construct each axis in place. We KNOW now that data_out has correct dimensions. We checked before getting here. We don't need to check len. It is the same as dims[i] each time. Perhaps we will anyway? :)
    tmp_axis = data_out->get_axis(i, len);
    if(len != dims[i]) return 1;

    N2 = ((float) dims[i])/2.0;
    res = this->get_res(i);
    for(int j= 0; j< dims[i]; j++) *(tmp_axis + j) = pi * ((float)j -N2)/N2/res;
  }

  return 0;

}

void data_array::copy_ids( data_array * src){
/** Copies ID fields from src array to this */

  strcpy(this->block_id, src->block_id);
  
  std::copy(src->time, src->time + 1, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src->space[i];
}

bool data_array::check_ids( data_array * src){
/** Checks ID fields match src */

  bool err=false;
  if(strcmp(this->block_id, src->block_id) != 0) err =true;
  
  for(int i=0; i< 3; i++) if(src->time[i] != this->time[i]) err=true;
  for(int i=0; i < 2; ++i) if(this->space[i] != src->space[i]) err=true;
  return err;
}


bool data_array::resize(int dim, int sz){
/** \brief Resize my_array on the fly
*
*dim is the dimension to resize, sz the new size. If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Similarly for axis elements. See my_array::resize() for more.
*/

  return 1;
  //call my_array::resize to resize data...
  bool err = my_array::resize(dim, sz);

  if(!err){
    //success! Do axes!
    my_type * new_ax;
    int part_sz = 0;

    if(dim == n_dims-1){
      //special case as we can shrink and maybe grow without copy
      for(int i=0; i<n_dims-1; ++i) part_sz+= dims[i];
      //product of all other dims
      new_ax = (my_type *) realloc((void*) this->axes, (part_sz+sz)*sizeof(my_type));
      if(!new_ax){
        my_print("Failed to reallocate axes", mpi_info.rank);
        return 1;
        //failure. leave as was.
      }
      int new_els = (sz - dims[n_dims-1]);
      
      if(new_els > 0) memset((void*)(new_ax + part_sz), 0.0, new_els*sizeof(my_type));
      //zero new elements
      axes = new_ax;
    }else{
      //have to allocate a new block and copy across.
      for(int i=0; i<dim-1; ++i) part_sz+= dims[i];
      for(int i=dim+1; i<n_dims; ++i) part_sz+= dims[i];

      new_ax=(my_type*)calloc(part_sz+sz,sizeof(my_type));
      int els_to_copy = 0, old_starts=0, new_starts=0;
      for(int i=0;i<n_dims; ++i){
        els_to_copy = dims[i];
        if(i == dim && sz< dims[i]) els_to_copy = sz;
        memcpy((void*)(axes + old_starts), (void*)(new_ax + new_starts), els_to_copy);

        old_starts += dims[i];
        i==dim ? new_starts += sz: new_starts += dims[i];
      
      }
      free(axes);
      axes = new_ax;

    }
    return 0;
  }else{
  
    return err;
  
  }


}


