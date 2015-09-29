//
//  main.cpp
//  
//
/** \file main.cpp \brief Main program
*
* This should: open a sequence of SDF files using the SDF library and read E and B field data. Fourier transform it. Extract frequency/wavenumber and angular spectra (if possible). Calculate the resulting particle diffusion coefficients using Lyons 1974 a, b, Albert 2005 and such.
  \author Heather Ratcliffe \date 18/09/2015.
*/
//
//

#include "main.h"
#include "support.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"

//#include <math>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>
//FFTW3 Fourier transform libs

using namespace std;

deck_constants my_const;

void abs_square( cplx_type * array, double * out, int nx);
void make_fft_axis(my_type * ax, int N, float res, int offset =0);
void test_bes();
spectrum * make_test_spectrum();

int main(int argc, char *argv[]){


cout<<"Code Version: "<<VERSION<<endl;

//TODO maybe we write a verify sdf which checks our files have the correct dimensionalities etc etc and contain needed blocks...

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
//TODO do something sensible like stop when the file doesn't exist and probably prompt to continue using data obtained...

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

//now we find the data for ex and make a data array with it, and axes
//for now, we assume we know it's 1-d so we make a 2-d array with 2 rows

int n_tims = 10;

data_array dat = data_array(block->dims[0], n_tims);
data_array dat_fft = data_array(block->dims[0], n_tims);

if(!dat.data or !dat_fft.data){
  cout<< "Bugger, data array allocation failed. Aborting."<<endl;
  return 0;
}
if(block->datatype_out != my_sdf_type){
  cout<< "Bugger wrong data type, recompile...";
  sdf_close(handle);
  return 0;
}

//dat.block_id = block->id;
//dat_fft.block_id = block->id;

strcpy(dat.block_id, block->id);
strcpy(dat_fft.block_id, block->id);
//set them to know what field they contain

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

int N = block->dims[0];

for(int i=0; i< n_tims; i++) dat.populate_row(block->data, block->dims[0], i);

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


//construct our FFT axes
//Hmm, lets put the original axis into our data array and the new data and axis into our result data array.

//Original axes are then copied from block->grids and constructed from the times. New are made by make_axis....

{
my_type * ax_ptr;
int len;

ax_ptr = dat.get_axis(0, len);
memcpy ((void *)ax_ptr, block->grids[0], len*sizeof(my_type));

dat.make_linear_axis(1, 1.0);
//generate uniform t axis of res 1
}

my_type * x_axis;
my_type * t_axis;

x_axis = (my_type*)malloc(N*sizeof(my_type));
t_axis = (my_type*)malloc(n_tims*sizeof(my_type));

make_fft_axis(x_axis, N, x_res);
make_fft_axis(t_axis, n_tims, 1.0);


sdf_close(handle);

//TODO this should all be wrapped away into a function which takes an input and output data_array and performs all this...

fftw_complex  *out;
double * in, *result;
float *result2;
fftw_plan p;


in = (double*) fftw_malloc(sizeof(double) * N*n_tims);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*n_tims);

std::copy(dat.data, dat.data+N*n_tims, in);
//copy because at the moment not a double... At some point we have to cast from the input type to double. May as well here.

p = fftw_plan_dft_r2c_2d(N, n_tims, in, out, FFTW_ESTIMATE);

fftw_execute(p);

result = (double*) fftw_malloc(sizeof(double) * N*n_tims);
result2 = (float*) fftw_malloc(sizeof(float) * N*n_tims);

abs_square(out, result, N*n_tims);


fftw_destroy_plan(p);
fftw_free(in);
fftw_free(out);

//test_bes();

std::copy(result, result + N*n_tims, result2);
//copy back to float. TODO fix this...

//copy into our data array, and add axes also
for(int i=0; i< n_tims; i++) dat_fft.populate_row(result2 + N*i, N, i);

{
my_type * ax_ptr;
int len;

ax_ptr = dat_fft.get_axis(0, len);
memcpy ((void *)ax_ptr, x_axis, len*sizeof(my_type));

ax_ptr = dat_fft.get_axis(1, len);
memcpy ((void *)ax_ptr, t_axis, len*sizeof(my_type));

}

free(x_axis);
free(t_axis);
fftw_free(result);
fftw_free(result2);

//Right now we have our FFT'd data with its bounds and axes in a decent structure.

//Next we write it to file to keep/visualise it

//bool my_array::write_to_file(std::fstream &file){
fstream file;
file.open("Tmp.txt", ios::out|ios::binary);

dat_fft.write_to_file(file);
file.close();

file.open("Tmp_tmp.txt", ios::in|ios::binary);

dat_fft.read_from_file(file);

file.close();



//Then we lineout and extract a wave spectrum
//Lets have a spectrum class then
//Contains data, axis, sizes, ids (field, time range, space range)
///HMMM. how to do angles?
//Relax seperability assumption so make it 2-d again....

//TODO modify spectrum to hold angle info.


spectrum * spect;
spect = make_test_spectrum();

//Then we use that and try and calculate the Diffusion coeff.




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

void make_fft_axis(my_type * ax, int N, float res, int offset){
//construct an fft'd axis of length N, from the original resolution. Units normed as input res.
//offset bevcause our axis array is one long consecutive 1-d one. Default is 0
//n_x2=float(n_pts)/2.
//ax=!pi*(findgen(n_pts)-n_x2)/n_x2/res

float N2;
N2 = ((float) N)/2.0;

//for(float i= -1* N2; i< N2; i++) *(ax + (int)i) = pi * i/N2/res;
for(int i= 0; i< N; i++) *(ax + i + offset) = pi * ((float)i - N2)/N2/res;


}

void test_bes(){
//test code using bessel funtion. Output values should be zeros.
//TODO expnd these things so we can test what we use before trying to proceed....


//cyl_bessel_j(v, x) = Jv(x)
//cyl_neumann(v, x) = Yv(x) = Nv(x)
cout<<"Bessel test: "<<endl;
double bess, arg;
int index;

index = 0;
arg = 2.40482555769577;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index = 1;
arg=7.01558666981561;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index=5;
arg=12.3386041974669;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess<<endl;

index =0;
arg =1.0;
bess = boost::math::cyl_bessel_j(index, arg);

cout<<bess-0.7651976865579665514497<<endl;

}

void get_deck_constants(){
//read deck.status and extract values for omega_ce, pe, vtherm etc
//TODO write this....

//NOTE any changes to the deck will need this modifying here....

my_const.omega_ce = 17588.200878;
my_const.omega_pe = 35176.401757;

}


spectrum * make_test_spectrum(){
//makes a basic spectrum object with suitable number of points, and twin, symmetric Gaussians centred at fixed x.

spectrum * ret;
char id[10];
id[0] = 'e'; id[1]='x';

ret = new spectrum(4096);
ret->set_ids(0, 100, 0, 4096, 1, id);

return ret;


}
