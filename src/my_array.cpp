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