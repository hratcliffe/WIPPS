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

#include "support.h"
#include "my_array.h"



my_array::my_array(int nx, int ny){
  defined = false;
  n_dims = 2;
  dims.push_back(nx);
  dims.push_back(ny);
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

bool my_array::populate_data(my_type * dat_in, int nx, int ny){
//needs to check type, check dimensions
//needs dimensions supplied...

return 1;
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
//TODO fill in in logical fashion....

file<<"Produced by "<<VERSION<<std::endl;
//dimension info
file<<n_dims<<" ";
for(int i=0;i<n_dims;i++) file<< dims[i]<<" ";
file<<std::endl;



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

