//
//  main.cpp
//  
//
//  Created by Heather Ratcliffe on 18/09/2015.
//
//

#include "main.h"
#include "my_array.h"
//#include <math>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <iostream>
#include "./SDF/C/include/sdf.h"
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]){

//Array class from std C++11 or boost
//own multi dim version
//version with axes


//Test use and syntax of FFTW

/*
my_array dat = my_array(5, 5);
bool err;
std::string ans;
ans = dat.array_self_test();
cout<<ans<<endl;
cout<<dat.get_element(3,6)<<endl;
cout<<dat.get_index(3,6)<<endl;
Tests my_array is working etc.
*/


int ierr, my_id, num_procs;
int err;
ierr = MPI_Init(&argc, &argv);
ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

sdf_file_t *handle = sdf_open("0000.sdf", MPI_COMM_WORLD, SDF_READ, 0);
//single threaded test!

if(handle){cout<<"Success!"<<endl;}
else{cout<<"Bleh"<<endl;}

err=sdf_read_blocklist(handle);
if(!err){cout<<"Yay!"<<endl;}
else{cout<<"Bum"<<endl;}

sdf_block_t * block, * next;

block = sdf_find_block_by_name(handle, "Electric Field/Ex");

if(block){cout<<"Success!!!"<<endl;}
else{cout<<"Bleh"<<endl;}

cout<<block->name<<" "<<block->id<<" "<<block->ndims<<" "<<block->dims[0]<<" "<<block->dims[1]<<" "<<block->dims[2]<<endl;


cout<<handle->time<<endl;
cout<<handle->nblocks<<endl;


next = handle->current_block;
for(int i =0; i< handle->nblocks; i++){
  cout<<next->name<<" "<<next->id<<endl;
  if(strcmp(next->name, "Electric Field/Ex")==0){
    cout<<"This one!"<<endl;
    cout<<"ndims "<<next->ndims<<endl;
    cout<<"dims are :"<<next->dims[0]<<" "<<next->dims[1]<<" "<<next->dims[2]<<endl;

  }

  //handle->current_block = next;
  next = next->next;
}

cout<<endl;
sdf_close(handle);


}
