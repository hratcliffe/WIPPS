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
#include "./SDF/C/include/stack_allocator.h"
#include "./SDF/C/include/sdf_helper.h"
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]){

//Array class from std C++11 or boost
//own multi dim version
//version with axes


//Test use and syntax of FFTW


//data_array dat = data_array(5, 5);
//cout<<"Good data array: "<<dat.is_good()<<endl;

/*
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

sdf_file_t *handle = sdf_open("0001.sdf", MPI_COMM_WORLD, SDF_READ, 0);
//single threaded test!

if(handle){cout<<"Success!"<<endl;}
else{cout<<"Bleh"<<endl;}

handle->stack_handle=stack_init();

err=sdf_read_blocklist(handle);
if(!err){cout<<"Yay!"<<endl;}
else{cout<<"Bum"<<endl;}

sdf_block_t * block, * next;

block = sdf_find_block_by_name(handle, "Electric Field/Ex");
block = sdf_find_block_by_id(handle, "ex");

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

//now we find the data for ex and make a data array with it, and axes
//for now, we assume we know it's 1-d so we make a 2-d array with 2 rows

data_array dat = data_array(block->dims[0], 2);

cout<<dat.get_dims(0)<<" "<<dat.get_dims(1)<<endl;

cout<<block->datatype_out<<" "<<SDF_DATATYPE_REAL4<<" "<<SDF_DATATYPE_REAL8<<endl;

cout<<block->dims[0];

if(block->data){cout<<"got data"<<endl;}
else{cout<<"uh OH!!!"<<endl;}
//now we memcpy from ex to our row
//but we need to match data type...

handle->current_block=block;
cout<<handle->current_block->name<<endl;

if(handle && block) cout<<"OK"<<endl;

handle->current_block = block;
//stack_alloc(handle->current_block);
sdf_read_data(handle);
//sdf_helper_read_data(handle, block);

cout<<"read"<<endl;
if(block->data){cout<<"got data"<<endl;}
else{cout<<"uh OH!!!"<<endl;}


if(block->datatype_out != my_sdf_type){
  cout<< "Bugger wrong data type, recompile...";
  sdf_close(handle);

  return 0;

}

float * my_ptr = (float *)block->data;

cout<<my_ptr[0]<<" "<<my_ptr[1]<<" "<<my_ptr[2]<<" "<<my_ptr[4095] <<endl;

//cout<<block->data[0]<<" "<<block->data[100]<<" "<<block->data[1023]<<" "<<block->data[4095]<<" "<<endl;

dat.populate_row(block->data, block->dims[0], 0);

cout<<dat.data[0]<<" "<<dat.data[1]<<" "<<dat.data[2]<<" "<<dat.data[4095]<<" "<<endl;

dat.populate_row(block->data, block->dims[0], 1);

cout<<dat.data[4096+0]<<" "<<dat.data[4096+1]<<" "<<dat.data[4096+2]<<" "<<dat.data[4096+4095]<<" "<<endl;

cout<<dat.get_element(0, 0)<<" "<<dat.get_element(0,1)<<endl;


sdf_close(handle);



}
