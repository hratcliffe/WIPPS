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
#include "support.h"
#include "my_array.h"
#include <complex.h>
#include <fftw3.h>
#include <cmath>


my_array::my_array(int nx, int ny){
  defined = false;
  n_dims = 2;
//  dims.push_back(nx);
//  dims.push_back(ny);
  dims = (int*)malloc(n_dims*sizeof(int));
  dims[0]=nx;
  dims[1]=ny;
  data = NULL;

  data=(my_type*)malloc(nx*ny*sizeof(my_type));

  if(data) defined = true;
  //Check allocation suceeded
}
my_array::~my_array(){

  if(data) free(data);
  data = NULL; // technically unnecessary as desctructor deletes memebers. 

}

my_type * my_array::get_ptr(int nx, int ny){

  int ind = get_index(nx, ny);
  if(ind  != -1){
    return data+ ind; //*sizeof(my_type);
  }else{
    return NULL;
  }

}

int my_array::get_index(int nx, int ny){
  
  if(n_dims != 2) return 1;
  if(nx < dims[0] && ny<dims[1]){
    return ny*dims[0] + nx;
  }else{
    return -1;
  }
}
int my_array::get_dims(){return n_dims;}
int my_array::get_dims(int dim){
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

my_type my_array::get_element(int nx, int ny){

  int ind = get_index(nx, ny);
  if(ind  != -1){
    return data[get_index(nx, ny)];
  }else{
    return 0;
  }

}

bool my_array::set_element(int nx, int ny, int val){
  //sets elements at nx, ny, and returns 1 if out of range, wrong no args, 0 else
  if(n_dims != 2) return 1;

  if(nx < dims[0] && ny<dims[1]){
    data[get_index(nx, ny)] = val;
    return 0;
  }else{
    return 1;
  }
  



}

bool my_array::populate_data(my_type * dat_in, int n_tot){
//Populates data with the total number of elements specified, as long as n_tot is less than the product of dims

int tot_els=1;
for(int i=0; i<n_dims;++i) tot_els *=dims[i];
if(n_tot > tot_els) return 1;

void * tmp = (void *) data;

if(!tmp) return 1;
memcpy (tmp , dat_in, n_tot*sizeof(my_type));

return 0;

}

bool my_array::populate_row(void * dat_in, int nx, int y_row){
//needs to check type, check dimensions
//needs dimensions supplied...

if(nx != dims[0]) return 1;
if(y_row > dims[1] || y_row < 0) return 1;

void * tmp = (void *)  get_ptr(0, y_row);

if(!tmp) return 1;
memcpy (tmp , dat_in, nx*sizeof(my_type));

return 0;


}

bool my_array::write_to_file(std::fstream &file){
//none of this will be human readable...
//Takes the version etc info then the whole data array and writes as a stream. It's not portable to other machines necessarily due to float sizes and endianness. It'll do. We'll start the file with a known float for confirmation.

/**IMPORTANT: this VERSION specifier links output files to code. If modifying output or order commit and clean build before using.
*/
const char tmp_vers[15] = VERSION;
char ch = '\0';

file.write((char*) &io_verify, sizeof(my_type));
file.write((char*) &tmp_vers, sizeof(char)*15);
//file.write(&ch, sizeof(char));

//Code version...

//dimension info
file.write((char*) &n_dims, sizeof(int));
int dim_tmp;
for(int i=0;i<n_dims;i++){
  dim_tmp = dims[i];
  file.write((char*) &dim_tmp, sizeof(int));
}
file.write((char *) data , sizeof(my_type)*dims[0]*dims[1]);

return 0;


}

bool my_array::read_from_file(std::fstream &file){
//For now this just spams to screen tested against what we wrote.

char tmp_vers[15];
my_type verf=0.0;

int n_dims_in, dim_tmp;

file.read((char*) &verf, sizeof(my_type));
std::cout<< verf<<" "<<io_verify<<std::endl;

file.read((char*) &tmp_vers, sizeof(char)*15);
std::cout<<tmp_vers<<" "<<VERSION<<std::endl;

if(verf != io_verify){
//equality even though floats as should be identical
  std::cout<<"Bugger, file read error";
  if(tmp_vers !=VERSION) std::cout<<"Incompatible code versions"<<std::endl;
  return 1;
}else{
  if(tmp_vers !=VERSION) std::cout<<"WARNING: A different code version was used to write this data. Proceed with caution. Fields may not align correctly."<<std::endl;
}


file.read((char*) &n_dims_in, sizeof(int));
std::cout<<n_dims_in<<" "<<n_dims<<std::endl;


for(int i=0;i<n_dims_in;i++){
  file.read((char*) &dim_tmp, sizeof(int));
  std::cout<<dim_tmp<<" "<<dims[i]<<std::endl;

}

my_type * data_tmp;
data_tmp=(my_type*)malloc(dims[0]*dims[1]*sizeof(my_type));

file.read((char *) data_tmp , sizeof(my_type)*dims[0]*dims[1]);
std::cout<<std::setprecision(9);
std::cout<<data_tmp[0]<<" "<<data[0]<<std::endl;
std::cout<<data_tmp[10]<<" "<<data[10]<<std::endl;
std::cout<<data_tmp[dims[0]-1]<<" "<<data[dims[0]-1]<<std::endl;
std::cout<<data_tmp[dims[0]+50]<<" "<<data[dims[0]+50]<<std::endl;
std::cout<<data_tmp[dims[0]*2+50]<<" "<<data[dims[0]*2+50]<<std::endl;

std::cout<<data_tmp[dims[0]*dims[1]-1]<<" "<<data[dims[0]*dims[1]-1]<<std::endl;

free(data_tmp);

return 0;


}


std::string my_array::array_self_test(){

bool err;
std::string ret;
int val;
//FOR 2-D only...
if(n_dims != 2) return "Only 2-D version exists";
//assign each element to unique val

for(int i=0; i<dims[0]; i++){
  for(int j =0; j<dims[1]; j++){
    err=set_element(i, j, (i+1)*(2*j+1));
    if(err) ret = "Cannot set element";
  }
}

//test assignments worked

for(int i=0; i<dims[0]; i++){
  for(int j =0; j<dims[1]; j++){
    val = get_element(i,j);
    if(val != (i+1)*(2*j+1)) ret += " Wrong element read";
    std::cout<<val<<" ";
  }
  std::cout<<std::endl;
}

if(ret == "") ret = "Array OK";

return ret;



}


data_array::data_array(int nx, int ny) : my_array(nx,ny){
//Constructor calls constructor for my_array and adds its own axes
  ax_defined = false;
  axes=(my_type*)malloc((nx+ny)*sizeof(my_type));
  if(axes) ax_defined=true;
}

data_array::~data_array(){
//Similarly destructor automatically calls destructor for my_array and frees axes

  if(axes) free(axes);
  axes = NULL; // technically unnecessary as desctructor deletes memebers.

}

my_type * data_array::get_axis(int dim, int & length){

  if(dim == 0){
    length = dims[dim];
    return axes;
  
  }else if(dim < n_dims){
    int offset=0;
    for(int i=0; i< dim; i++) offset +=dims[i];
    length = dims[dim];
    return axes + offset;
  
  }else{
    length = 0;
    return NULL;
  }

}

void data_array::make_linear_axis(int dim, float res, int offset){

int len;
my_type * ax_ptr = get_axis(dim, len);

for(int i=0; i<len; i++){
  *(ax_ptr +i) = ((float) (i-offset)) * res;

}


}

bool data_array::write_to_file(std::fstream &file){

/**IMPORTANT: the VERSION specifier links output files to code. If modifying output or order commit and clean build before using.
*/

char ch = '\0';
file.write(block_id, sizeof(char)*10);
file.write(&ch, sizeof(char));
//Null terminating character

my_array::write_to_file(file);
//call base class method to write that data.

file.write((char *) axes , sizeof(my_type)*(dims[0]+dims[1]));
//Add axes.

return 0;

}

bool data_array::read_from_file(std::fstream &file){
//For now this just spams to screen tested against what we wrote.

//First read the block ID
char id_in[11];

file.read(id_in, sizeof(char)*11);
std::cout<< id_in<<" "<<block_id<<std::endl;

my_array::read_from_file(file);
//call parent class to read data

//now read axes
std::cout<<"Axes :"<<std::endl;
my_type * data_tmp;
data_tmp=(my_type*)malloc((dims[0]+dims[1])*sizeof(my_type));

file.read((char *) data_tmp , sizeof(my_type)*(dims[0]+dims[1]));

std::cout<<data_tmp[0]<<" "<<axes[0]<<std::endl;
std::cout<<data_tmp[10]<<" "<<axes[10]<<std::endl;
std::cout<<data_tmp[dims[0]-1]<<" "<<axes[dims[0]-1]<<std::endl;
std::cout<<data_tmp[dims[0]]<<" "<<axes[dims[0]]<<std::endl;
std::cout<<data_tmp[dims[0]+dims[1]-1]<<" "<<axes[dims[0]+dims[1]-1]<<std::endl;

free(data_tmp);


return 0;


}

bool data_array::fft_me(data_array * data_out){
/** \brief FFT data_array
*
* Data and axes in this object are FFT'd using FFTW and stored into the instance pointed to by data_out. Data_out must be created with correct dimensions first, but we check and return error (1) if it is not so.
*/

if(!data_out->is_good()){
  std::cout<<"Output array for FFT undefined"<<std::endl;
  return 1;
}
if(data_out->n_dims != this->n_dims){
  std::cout<<"Wrong output dimensions for FFT"<<std::endl;
  return 1;
}
for(int i=0; i<n_dims;++i){
  if(data_out->dims[i] != this->dims[i]){
    std::cout<<"Wrong output dimensions for FFT"<<std::endl;
    return 1;
  }
}

int total_size=1; /**< Total number of elements in array*/
for(int i=0; i<n_dims;++i) total_size *= dims[i];

int fft_dim =1;/**< Dimension to FFT over, if required*/

ADD_FFTW(plan) p;
cplx_type *out;
my_type * in, *result;

in = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);
//my_type should match the used FFTW library, so no type conversion necessary
out = (cplx_type *) fftwf_malloc(sizeof(cplx_type) * total_size);

result = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);

//Possibly this bit can be genericised?
if(n_dims == 1){
  p = ADD_FFTW(plan_dft_r2c_1d)(dims[0], in, out, FFTW_ESTIMATE);

}else if(n_dims == 2){
  p = ADD_FFTW(plan_dft_r2c_2d)(dims[0], dims[1], in, out, FFTW_ESTIMATE);

}else{
  return 1;
}

ADD_FFTW(execute)(p);
//Execute the plan

cplx_type * addr;
addr = out;
//because double indirection is messy and cplx type is currently a 2-element array of floats

for(int i=0; i< total_size; i++){
  *(result+i) = (my_type)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1])) ;
  addr++;
}
//Absolute square of out array to produce final result of type my_type

bool err;
err = data_out->populate_data(result, total_size);
//Copy result into out array
std::cout<<err<<std::endl;

ADD_FFTW(destroy_plan)(p);
ADD_FFTW(free)(in);
ADD_FFTW(free)(out);
ADD_FFTW(free)(result);
//Destroy stuff we don't need


my_type * tmp_axis;
float N2, res;
int len;

for(int i=0;i<n_dims;++i){
//Loop over the dimensions and construct each axis in place. We KNOW now that data_out has correct dimensions. We checked before getting here. We don't need to check len. It is the same as dims[i] each time. We will anyway :)
  N2 = ((float) dims[i])/2.0;
  res = this->get_res(i);
  tmp_axis = data_out->get_axis(i, len);
  if(len != dims[i]) return 1;
  for(int j= 0; j< dims[i]; i++) *(tmp_axis + j) = pi * ((float)j - N2)/N2/res;
}

return 0;

}


float data_array::get_res(int i){
//return resolution of axis on dimension i. Assumes linear etc etc
int len;
my_type * axis = this->get_axis(i, len);

return std::abs(axis[0]-axis[1]);


}

