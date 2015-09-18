//
//  main.h
//  
//
//  Created by Heather Ratcliffe on 18/09/2015.
//
//

#ifndef ____main__
#define ____main__

#include <stdio.h>
#include <vector>
#include <fstream>

#define my_type float
//Shortcut so we can change array internal type later

using namespace std;
class my_array{

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

}

bool is_good(){return defined;}

my_type get_element(int nx, int ny){
  return data[get_index(nx, ny)];

}
int get_index(int nx, int ny){

  return nx*dims[1] + nx;

}


};


#endif /* defined(____main__) */
