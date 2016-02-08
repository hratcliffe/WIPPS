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
#include <cmath>
#include "stripped_array.h"


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
  n_dims = 0;
  data=nullptr;

}

my_array::my_array(long nx, long ny){
/** \brief 2-d rectangular array
*
*Sets up a n-d rectangular array of nx*ny and allocates data arrays. If either of first 2 sizes is zero or exceeds MAX_SIZE, print error and exit. Otherwise construct 2, 3, 4d as dims specify
*/
  construct();

  if(ny==0 || nx==0){
    return;
  }
  if(ny> MAX_SIZE || nx>MAX_SIZE){
    return;
  }
  std::cout<<"Constructifying"<<std::endl;

  n_dims = 2;
  this->dims = (long*)malloc(n_dims*sizeof(long));
  dims[0]=nx;
  dims[1]=ny;

  data=(my_type*)calloc(nx*ny,sizeof(my_type));
  
  if(data){
    defined = true;
  }else{
  
    std::cout<<"Cannot allocate"<<std::endl;
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
}

long my_array::get_index(long nx, long ny){

   // if((nx < dims[0]) && (ny<dims[1])){
  // int ret;
      return (((nx < dims[0]) && (ny<dims[1])) ?  ny*dims[0] + nx : -1);
//      return ret;
    //}else{
      //return -1;
   // }
}

int my_array::get_dims(){
  return n_dims;
}
int my_array::get_dims(int dim){
  if(dim < n_dims){
    return dims[dim];
  
  }else{return 0;}
}

my_type my_array::get_element(int nx, int ny){
/** Return element at nx, ny. Out of range etc will return 0.0*/
  int ind = get_index(nx, ny);

/*  if(ind  != -1){
    return data[ind];
  }else{
    return 0.0;
  }*/
  return ((ind==-1) ? 0.0 : data[ind] );

}

bool my_array::set_element(long nx, long ny, my_type val){
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
