//
//  my_array.h
//  
//
//  Created by Heather Ratcliffe on 21/09/2015.
//
//

#ifndef ____my_array__
#define ____my_array__

#include <stdio.h>
#include <stdlib.h> //for malloc eeejit
#include <vector>

#define my_type float
//Shortcut so we can change array internal type later if needed

using namespace std;
class my_array{
protected:
int n_dims;
std::vector<int> dims;
bool defined;
//flag to check memory allocation sucess

public:

my_type *data;

my_array(int nx, int ny){
  defined = false;
  n_dims = 2;
  dims.push_back(nx);
  dims.push_back(ny);
  data = NULL;

  data=(my_type*)malloc(nx*ny*sizeof(my_type));

  if(data) defined = true;
  //Check allocation suceeded
}

bool is_good(){return defined;}

my_type get_element(int nx, int ny){
  return data[get_index(nx, ny)];

}
int get_index(int nx, int ny){

  return nx*dims[1] + nx;

}
int get_dims(){return n_dims;}

//we can use [] to wrap get elements and have the args pushed into vector which we then work with to be generic

};



class data_array : protected my_array{
//Adds axes for data

bool ax_defined;
public:

my_type *axes;
//This will be 1-d array in sections, so can be arbitary length and dims


data_array(int nx, int ny) : my_array(nx,ny){
  ax_defined = false;
  axes=(my_type*)malloc((nx+ny)*sizeof(my_type));
  if(axes) ax_defined=true;
}

bool is_good(){return defined && ax_defined;}

my_type * get_axis(int dim, int & length){

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

};



#endif /* defined(____my_array__) */
