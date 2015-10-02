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

int n_tims = max(tim_in[1]-tim_in[0], 1);

int n_dims;
std::vector<int> dims;
my_reader.read_dims(n_dims, dims);

if(n_dims !=1) return 1;
//for now abort if data file wrong size...

data_array dat = data_array(dims[0], n_tims);
data_array dat_fft = data_array(dims[0], n_tims);

if(!dat.data or !dat_fft.data){
  cout<< "Bugger, data array allocation failed. Aborting."<<endl;
  return 0;
}
//Make arrays and check success.

my_reader.read_data(&dat, tim_in, space_in);
//read raw data and axes

/*
float * tmp_dat;

tmp_dat = (float*)malloc(N*sizeof(float));
for(int i=0; i<N ; i++) *(tmp_dat+i) = sin((float)i /25.0);
//sin(x/2500.0)
//Now we generate some sine data instead to test FFT

for(int i=0; i< n_tims; i++) dat.populate_row(tmp_dat, N, i);

*/

err = dat.fft_me(&dat_fft);

cout<<"FFT returned err_state "<<err<<endl;


//Right now we have our FFT'd data with its bounds and axes in a decent structure.
//Next we write it to file to keep/visualise it

fstream file;
file.open("Tmp.txt", ios::out|ios::binary);

dat.write_to_file(file);
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
