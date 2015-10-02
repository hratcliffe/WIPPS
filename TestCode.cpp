
//Contains all the temporary and test code fragments that I might want later.

//data_array dat = data_array(5, 5);
//cout<<"Good data array: "<<dat.is_good()<<endl;

/*
ans = dat.array_self_test();
cout<<ans<<endl;
cout<<dat.get_element(3,6)<<endl;
cout<<dat.get_index(3,6)<<endl;
Tests my_array is working etc.
*/


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


//These bits print some stuff so we check the copying worked...
/*cout<<my_ptr[0]<<" "<<my_ptr[1]<<" "<<my_ptr[2]<<" "<<my_ptr[4095] <<endl;

cout<<dat.data[0]<<" "<<dat.data[1]<<" "<<dat.data[2]<<" "<<dat.data[4095]<<" "<<endl;

cout<<dat.data[4096+0]<<" "<<dat.data[4096+1]<<" "<<dat.data[4096+2]<<" "<<dat.data[4096+4095]<<" "<<endl;

cout<<dat.get_element(0, 0)<<" "<<dat.get_element(0,1)<<endl;
*/

/* 
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
*/

/*
#include <stdint.h>
#if UINTPTR_MAX == 0xffffffff
/ 32-bit /
#elif UINTPTR_MAX == 0xffffffffffffffff
/ 64-bit /
#else
/ wtf /
#endif

*/


//Test code to read back vals and check against originals
/*file.open("Tmp.txt", ios::in|ios::binary);
my_type tmp, tmp2;

for(int i=0; i<4096; i++){
  tmp2 = dat_fft.get_element(i,0);
  file.read((char *) &tmp, sizeof(my_type));
  cout<<tmp-tmp2<<" ";

}*/

/*file<<n_dims<<" ";
for(int i=0;i<n_dims;i++) file<< dims[i]<<" ";
file<<std::endl;

file.write((char *) data , sizeof(my_type)*dims[0]);
*/

/*
fftwf_complex  *out;
float * in, *result;
my_type *result2;
fftwf_plan p;

int dims[4];

dims[0] = N;
dims[1] = n_tims;
//in = (double*) fftw_malloc(sizeof(double) * N*n_tims);
in = (float*) fftwf_malloc(sizeof(float) * N*n_tims);
out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N*n_tims);

//copy because at the moment not a double... At some point we have to cast from the input type to double. May as well here.
//Check worked OK

cout<<"Checking type promotion"<<endl;
cout<< in[0]<<" "<<in[N]<<" "<<dat.get_element(0,0)<<" "<<dat.get_element(0,1)<<endl;
cout<< in[50]<<" "<<in[N+50]<<" "<<dat.get_element(50,0)<<" "<<dat.get_element(50,1)<<endl;


//p = fftw_plan_dft_r2c_2d(N, n_tims, in, out, FFTW_ESTIMATE);
//p = fftwf_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

int n_dim =1;
//Dimension for FFT
p = fftwf_plan_dft_r2c(n_dim, dims, in,out, FFTW_MEASURE);

std::copy(dat.data, dat.data+N*n_tims, in);

//fftwf_execute_dft_r2c(p, in, out);
fftwf_execute(p);

//result = (double*) fftw_malloc(sizeof(double) * N*n_tims);
result2 = (my_type*) fftwf_malloc(sizeof(my_type) * N*n_tims);
result = (float*) fftwf_malloc(sizeof(float) * N*n_tims);


abs_square(out, result2, N*n_tims);

fftwf_destroy_plan(p);
fftwf_free(in);
fftwf_free(out);

//test_bes();

//std::copy(result, result + N*n_tims, result2);
//copy back to float. TODO fix this...

cout<<"And demotion"<<endl;
cout<< result[0]<<" "<<result2[0]<<endl;
cout<<result[N]<<" "<<result2[N]<<endl;
cout<< result[50]<<" "<<result2[50]<<endl;
cout<<result[N+50]<<result2[N+50]<<" "<<endl;


//copy into our data array, and add axes also
for(int i=0; i< n_tims; i++) dat_fft.populate_row(result2 + N*i, N, i);

dat_fft.populate_row(result2, N, 0);

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
fftwf_free(result);
fftwf_free(result2);
*/

