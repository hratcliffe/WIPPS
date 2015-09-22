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
#define my_sdf_type 3
//Shortcut so we can change array internal type later if needed. And matching SDF type

template<class T> class arrayl{

T* array_values;
T getvalue(int ix,int iy);

};

//using namespace std;
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

~my_array(){

  if(data) free(data);
  data = NULL; // technically unnecessary as desctructor deletes memebers. 

}

bool is_good(){return !defined;}

my_type get_element(int nx, int ny){

  int ind = get_index(nx, ny);
  if(ind  != -1){
    return data[get_index(nx, ny)];
  }else{
    return 0;
  }

}

my_type * get_ptr(int nx, int ny){

  int ind = get_index(nx, ny);
  if(ind  != -1){
    return data; //+ ind*sizeof(my_type);
  }else{
    return 0;
  }

}

int get_index(int nx, int ny){
  
  if(n_dims != 2) return 1;
  if(nx < dims[0] && ny<dims[1]){
    return ny*dims[0] + nx;
  }else{
    return -1;
  }
}
int get_dims(){return n_dims;}
int get_dims(int dim){
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

bool set_element(int nx, int ny, int val){
  //sets elements at nx, ny, and returns 1 if out of range, wrong no args, 0 else
  if(n_dims != 2) return 1;

  if(nx < dims[0] && ny<dims[1]){
    data[get_index(nx, ny)] = val;
    return 0;
  }else{
    return 1;
  }
  



}

bool populate_data(my_type * dat_in, int nx, int ny);
bool populate_row(void * dat_in, int nx, int y_row);


//we can use [] to wrap get elements and have the args pushed into vector which we then work with to be generic

std::string array_self_test();

};



class data_array : public my_array{
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
