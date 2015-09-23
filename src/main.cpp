//
//  main.cpp
//  
//
//  Created by Heather Ratcliffe on 18/09/2015.
//
//

#include "../include/main.h"
#include "../include/my_array.h"
//#include <math>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <iostream>
#include "../SDF/C/include/sdf.h"
#include "../SDF/C/include/stack_allocator.h"
#include "../SDF/C/include/sdf_helper.h"
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>

using namespace std;

#define cplx_type fftw_complex

const float pi = 3.14159;

void abs_square( cplx_type * array, double * out, int nx);
void make_axis(my_type * ax, int N, float res, int offset =0);

int main(int argc, char *argv[]){

//Array class from std C++11 or boost
//own multi dim version
//version with axes


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

if(handle){cout<<"Success! File opened"<<endl;}
else{
  cout<<"Bleh. File open error. Aborting"<<endl;
  return 0;
}

//handle->stack_handle=stack_init();

err=sdf_read_blocklist(handle);
if(!err){cout<<"Yay! Blocks read"<<endl;}
else{cout<<"Bum. Block read failed"<<endl;}

sdf_block_t * block, * next;

//block = sdf_find_block_by_name(handle, "Electric Field/Ex");
block = sdf_find_block_by_id(handle, "ex");

if(block){cout<<"Success!!! Ex found"<<endl;}
else{
  cout<<"Bleh. No Ex found. Aborting"<<endl;
  return 0;
}

cout<<block->name<<" "<<block->id<<" "<<block->ndims<<" "<<block->dims[0]<<" "<<block->dims[1]<<" "<<block->dims[2]<<endl;

cout<<"File time "<<handle->time<<endl;
//cout<<handle->nblocks<<endl;

/*
//Cycles though blocks and prints their names, and extra info for the ex one
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
*/

//now we find the data for ex and make a data array with it, and axes
//for now, we assume we know it's 1-d so we make a 2-d array with 2 rows

data_array dat = data_array(block->dims[0], 2);

if(!dat.data){
  cout<< "Bugger, data array allocation failed. Aborting."<<endl;
  return 0;
}
if(block->datatype_out != my_sdf_type){
  cout<< "Bugger wrong data type, recompile...";
  sdf_close(handle);
  return 0;
}

//cout<<block->datatype_out<<" "<<SDF_DATATYPE_REAL4<<" "<<SDF_DATATYPE_REAL8<<endl;

handle->current_block = block;
//stack_alloc(handle->current_block);
sdf_read_data(handle);
//sdf_helper_read_data(handle, block);

if(block->data){cout<<"Got data"<<endl;}
else{
  cout<<"uh OH!!! Read failed. Aborting"<<endl;
  return 0;
}


float * my_ptr = (float *)block->data;

dat.populate_row(block->data, block->dims[0], 0);
dat.populate_row(block->data, block->dims[0], 1);

//These bits print some stuff so we check the copying worked...
cout<<my_ptr[0]<<" "<<my_ptr[1]<<" "<<my_ptr[2]<<" "<<my_ptr[4095] <<endl;

cout<<dat.data[0]<<" "<<dat.data[1]<<" "<<dat.data[2]<<" "<<dat.data[4095]<<" "<<endl;

cout<<dat.data[4096+0]<<" "<<dat.data[4096+1]<<" "<<dat.data[4096+2]<<" "<<dat.data[4096+4095]<<" "<<endl;

cout<<dat.get_element(0, 0)<<" "<<dat.get_element(0,1)<<endl;

//now get the grids res.
//and also grab two consecutive times to check t res?
float x_res=0;

block = sdf_find_block_by_id(handle, "grid");

//checks type is plain i/e/ even gridded mesh
if(block->blocktype != SDF_BLOCKTYPE_PLAIN_MESH){cout<< "Uh oh, grids look wrong"<<endl;}

handle->current_block = block;
sdf_read_data(handle);
my_ptr = (float *)block->grids[0];
if(my_ptr){x_res = my_ptr[1] - my_ptr[0];}
cout<<"x resolution is "<<x_res<<endl;
//gets x-axis resolution


sdf_close(handle);

//temproary hard definitions
int N = 4096;
fftw_complex  *out;
double * in, *result;
fftw_plan p;

in = (double*) fftw_malloc(sizeof(double) * N);

std::copy(dat.data, dat.data+4096, in);

out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

p = fftw_plan_dft_r2c_1d(N, in, out,FFTW_ESTIMATE);
//FFTW plans find best way to perform the FFT before doing it. Ideal for when have multiple ones to do. Estimate is less optimised, but quicker to find. If have many to do, can use FFTW_CALCULATE
//fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

fftw_execute(p);

result = (double*) fftw_malloc(sizeof(double) * N);
abs_square(out, result, N);
//copies out to result also

fftw_destroy_plan(p);
fftw_free(in);
fftw_free(out);


//result is now an FFT'd Ex array
//construct our axes
my_type * x_axis;

x_axis = (my_type*)malloc(N*sizeof(my_type));

//got x grid res from sdf block above

make_axis(x_axis, N, x_res);


free(x_axis);
fftw_free(result);

}

//elements wise ops so treat as long 1-d arrray. nx should be total length product(dims). Input is pointer, needs to be to pre-defined memory of apporpireate size. We'll move all this inside our arrays later
void abs_square( cplx_type * array, double * out, int nx){

cplx_type * addr;
//because double indirection is messy and cplx type is currently a 2-element array of doubles

for(int i=0; i< nx; i++){
  addr = (array + i);
  *(out+i) = (*addr[0])*(*addr[0]) + (*addr[1])*(*addr[1]) ;
}

}

void make_axis(my_type * ax, int N, float res, int offset){
//construct an fft'd axis of length N, from the original resolution. Units normed as input res.
//offset bevcause our axis array is one long consecutive 1-d one. Default is 0
//n_x2=float(n_pts)/2.
//ax=!pi*(findgen(n_pts)-n_x2)/n_x2/res

float N2;
N2 = ((float) N)/2.0;

//for(float i= -1* N2; i< N2; i++) *(ax + (int)i) = pi * i/N2/res;
for(int i= 0; i< N; i++) *(ax + i + offset) = pi * ((float)i - N2)/N2/res;


}

