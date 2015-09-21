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

#include "my_array.h"



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