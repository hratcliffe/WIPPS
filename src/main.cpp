/** \file main.cpp \brief Main program
*
* This should: open a sequence of SDF files using the SDF library and read E and B field data. Fourier transform it. Extract frequency/wavenumber and angular spectra (if possible). Calculate the resulting particle diffusion coefficients using Lyons 1974 a, b, Albert 2005 and such.
  \author Heather Ratcliffe \date 18/09/2015.
*/

#include "main.h"
#include "support.h"
#include "reader.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"

#include <math.h>
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

deck_constants my_const;/*< Physical constants*/


void abs_square( cplx_type * array, double * out, int nx);
void abs_square( cplx_type * array, float * out, int nx);
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


char block_id[10] = "ex";
reader my_reader("", block_id);
//reader(std::string file_prefix_in,  char * block_id_in){

int tim_in[2], space_in[2];
tim_in[0]=0;
tim_in[1]=10;
space_in[0]=0;
space_in[1]=-1;



sdf_file_t *handle = sdf_open("0001.sdf", MPI_COMM_WORLD, SDF_READ, 0);
//single threaded test!

if(handle){cout<<"Success! File opened"<<endl;}
else{
  cout<<"Bleh. File open error. Aborting"<<endl;
  return 0;
}
//TODO do something sensible like stop when the file doesn't exist and probably prompt to continue using data obtained...

err=sdf_read_blocklist(handle);
if(!err){cout<<"Yay! Blocks read"<<endl;}
else{cout<<"Bum. Block read failed"<<endl;}

sdf_block_t * block, * next;

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
//for now, we assume we know it's 1-d so we make a 2-d array with duplicated rows

int n_tims = 1;

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

strcpy(dat.block_id, block->id);
strcpy(dat_fft.block_id, block->id);
//set them to know what field they contain


my_reader.read_data(&dat, tim_in, space_in);
//bool read_data(data_array * my_data_in, int time_range[2], int space_range[2]);


handle->current_block = block;
sdf_read_data(handle);

//sdf_close(handle);


if(block->data){cout<<"Got data"<<endl;}
else{
  cout<<"uh OH!!! Read failed. Aborting"<<endl;
  return 0;
}


float * my_ptr = (float *)block->data;

int N = block->dims[0];

//TODO this only works if my_type is a float...
bool err2 =0;
for(int i=0; i< n_tims; i++){

  err2 = dat.populate_row(block->data, N, i);
  if(err2) cout<<"Bugger!"<<endl;
}

//Check data in row matches up...
cout<<"Chekcing"<<endl;
cout<< my_ptr[0]<<" "<<dat.get_element(0,0)<<" "<<dat.get_element(0,1)<<endl;
cout<< my_ptr[50]<<" "<<dat.get_element(50,0)<<" "<<dat.get_element(50,1)<<endl;

float * tmp_dat;

tmp_dat = (float*)malloc(N*sizeof(float));
for(int i=0; i<N ; i++) *(tmp_dat+i) = sin((float)i /25.0);
//sin(x/2500.0)
//Now we generate some sine data instead to test FFT

for(int i=0; i< n_tims; i++) dat.populate_row(tmp_dat, N, i);

cout<<" "<<dat.get_element(0,0)<<" "<<dat.get_element(100,0)<<" "<<dat.get_element(2048,0)<<" "<<dat.get_element(4095,0)<<endl;


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


{
my_type * ax_ptr;
int len;

ax_ptr = dat.get_axis(0, len);
//TODO this only works right is my_type is float...
memcpy ((void *)ax_ptr, block->grids[0], len*sizeof(my_type));

dat.make_linear_axis(1, 1.0);
//generate uniform t axis of res 1
}


sdf_close(handle);

err = dat.fft_me(&dat_fft);

cout<<"FFT returned err_state "<<err<<endl;


//Right now we have our FFT'd data with its bounds and axes in a decent structure.
//Next we write it to file to keep/visualise it

fstream file;
file.open("Tmp.txt", ios::out|ios::binary);

dat_fft.write_to_file(file);
//dat.write_to_file(file);

file.close();

//file.open("Tmp_dat.txt", ios::in|ios::binary);

//dat.read_from_file(file);

//file.close();



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

cplx_type * addr = array;
//because double indirection is messy and cplx type is currently a 2-element array of doubles (maybe floats. So cast explicitly)

for(int i=0; i< nx; i++){
//  addr = (array + i);
  *(out+i) = (double)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1])) ;

  addr++;

}

}

void abs_square( cplx_type * array, float * out, int nx){

cplx_type * addr = array;
//because double indirection is messy and cplx type is currently a 2-element array of floats

for(int i=0; i< nx; i++){
  *(out+i) = (float)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1])) ;
  addr++;

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
