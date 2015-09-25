


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