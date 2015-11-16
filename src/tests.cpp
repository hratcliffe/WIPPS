//
//  tests.cpp
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

/** \file tests.cpp 
*
*To add a test, do the following:
*Descend an object from test_entity which has at least a constructor doing any setup required, a name string for output id, a function run taking no parameters which performs the necessary test and a destructor doing cleanup. Add any other member variables or functions required. In tests::setup_tests create an instance of your class as test_obj = new your_class() and then add your test to the remit using add_test(test_obj); Alternately make the instance and use the global test_bed using test_bed->add(your pntr) from anywhere.
*To add errors, add the message into the blank spaces in the list below, err_names, and declare a const int of desired name aliased to TEST_USERDEF_ERR* where * = 1-4

*/

#include <stdio.h>
#include <math.h>
#include "tests.h"
#include "reader.h"
#include "support.h"
#include "my_array.h"

extern mpi_info_struc mpi_info;
extern tests * test_bed;
class reader;

const int err_codes[err_tot] ={TEST_PASSED, TEST_WRONG_RESULT, TEST_NULL_RESULT, TEST_ASSERT_FAIL, TEST_USERDEF_ERR1, TEST_USERDEF_ERR2, TEST_USERDEF_ERR3, TEST_USERDEF_ERR4};

std::string err_names[err_tot]={"None", "Wrong result", "Invalid Null result", "Assignment or assertion failed", "", "", "", ""};

tests::tests(){
  setup_tests();
}
tests::~tests(){
  cleanup_tests();
}

void tests::setup_tests(){
/** \brief Setup test bed
*
*Opens reporting file. Then instantiates all the test objects and adds them into the test_list
*/
  outfile = new std::fstream();
  outfile->open(filename.c_str(), std::ios::out);
  if(!outfile->is_open()){
    std::cout<<"Error opening "<<filename<<std::endl;
    //can't log so return with empty test list
    return;
  }

  test_entity * test_obj;

  test_obj = new test_entity_reader();
  add_test(test_obj);
  //these two lines are needed for each test you want to do.
  test_obj = new test_entity_data_array();
  add_test(test_obj);
  test_obj = new test_entity_get_and_fft();
  add_test(test_obj);
  test_obj = new test_entity_basic_maths();
  add_test(test_obj);

}
void tests::add_test(test_entity * test){
  /** Adds a test to the list to be performed*/
  test_list.push_back(test);
}

void tests::report_err(int err, int test_id){
  if(test_id == -1) test_id = current_test_id;

  my_print(outfile, get_printable_error(err, test_id), mpi_info.rank);

}

std::string tests::get_printable_error(int err, int test_id){

  std::string err_string="";
  if(err!=TEST_PASSED){
    for(int i=err_tot-1; i>0; --i){
      //Run most to least significant
      if((err & err_codes[i]) == err_codes[i]){
        err_string +=err_names[i] + ", ";
      }
    }
  
    err_string = "Error "+err_string+"(code "+mk_str(err)+") on";
  }
  else err_string = "Passed";
  return err_string+" test "+test_list[test_id]->name;

}

void tests::cleanup_tests(){
/** \brief Delete test objects
*
*
*/
  if(outfile->is_open()){
    std::cout<<"Testing complete and logged in "<<filename<<std::endl;
    outfile->close();
  }else{
    std::cout<<"No log generated"<<std::endl;
  }
  delete outfile;
  for(current_test_id=0; current_test_id<test_list.size(); current_test_id++){
    delete test_list[current_test_id];
  
  }

  
}

void tests::run_tests(){
/** \brief Delete test objects
*
*
*/
  int err = TEST_PASSED;
  for(current_test_id=0; current_test_id<test_list.size(); current_test_id++){
    test_list[current_test_id]->run();
    
  }

}

test_entity_reader::test_entity_reader(){
  name = "reader class";
  char block_id[10]= "run_info";
  test_rdr = new reader("test", block_id);
  
}
test_entity_reader::~test_entity_reader(){

  delete test_rdr;
}

int test_entity_reader::run(){
  int err = TEST_PASSED;
  
  int sz = test_rdr->get_file_size();
  if(sz != size) err |= TEST_WRONG_RESULT;

/* do testing */
  test_bed->report_err(err);
  return err;

}

test_entity_data_array::test_entity_data_array(){
  name = "data array class";
  test_array = new data_array(10, 10);
  
}
test_entity_data_array::~test_entity_data_array(){

  delete test_array;
}

int test_entity_data_array::run(){
  int err = TEST_PASSED;
  bool tmp_err;
  
  int val;
  if(test_array->get_dims() == 2){
    //assign each element to unique val

    for(int i=0; i<test_array->get_dims(0); i++){
      for(int j =0; j<test_array->get_dims(1); j++){
        tmp_err=test_array->set_element(i, j, (i+1)*(2*j+1));
        if(tmp_err) err |= TEST_ASSERT_FAIL;
      }
    }

    //test assignments worked

    for(int i=0; i<test_array->get_dims(0); i++){
      for(int j =0; j<test_array->get_dims(1); j++){
        val = test_array->get_element(i,j);
        if(val != (i+1)*(2*j+1)) err |=TEST_WRONG_RESULT;
      }
    }
  }

/* do testing */
  test_bed->report_err(err);
  return err;

}

test_entity_get_and_fft::test_entity_get_and_fft(){
  name = "read and FFT";
  char block_id[10]= "ex";
  test_rdr = new reader("test", block_id);

}
test_entity_get_and_fft::~test_entity_get_and_fft(){
  delete test_rdr;

}

int test_entity_get_and_fft::run(){

  int err = TEST_PASSED;

  int tim_in[2], space_in[2];
  tim_in[0]=0;
  tim_in[1]=1;
  space_in[0]=0;
  space_in[1]=-1;

  int n_tims = std::max(tim_in[1]-tim_in[0], 1);

  int n_dims;
  std::vector<int> dims;
  test_rdr->read_dims(n_dims, dims);
  std::cout<<dims[0]<<std::endl;
  if(n_dims !=1){
    err |= TEST_WRONG_RESULT;
    return err;
    //nothing more worth doing right now...
  }

  test_dat = new data_array(dims[0], n_tims);
  test_dat_fft = new data_array(dims[0], n_tims);
  if(!test_dat->data||!test_dat_fft->data){
    err|=TEST_ASSERT_FAIL;
    return err;
  }

  test_rdr->read_data(test_dat, tim_in, space_in);

  bool tmp_err = test_dat->fft_me(test_dat_fft);
  if(tmp_err) err|=TEST_ASSERT_FAIL;
  
  /** Now test the resulting frequency is right.... Also tests our axes...*/
  
  //Get primary frequency
  
  test_bed->report_err(err);
  return err;

}

test_entity_basic_maths::test_entity_basic_maths(){

  name = "basic maths helpers";
  size = 500;
  data_square=(calc_type*)calloc(size,sizeof(calc_type));
  data_positive=(calc_type*)calloc(size,sizeof(calc_type));
  data_tmp=(calc_type*)calloc(size,sizeof(calc_type));
  axis=(calc_type*)calloc(size,sizeof(calc_type));
  d_axis=(calc_type*)calloc(size,sizeof(calc_type));

  data_square[0] = 1.0;
  data_positive[0] = 0.0;
  axis[0] = 0.0;
  d_axis[0] = 1.0;
  for(int i=1; i<size; ++i){
    data_square[i] = - data_square[i-1];
    //alternating square wave, average=int=0
    data_positive[i] = data_positive[i-1] + 0.1;
    //monotonic increase, integral = size * size/10/2 depending on upper bnd
    d_axis[i] = 1.0;
    axis[i] = axis[i-1] + d_axis[i];
  }

//set up some data arrays...

}
test_entity_basic_maths::~test_entity_basic_maths(){

  free(data_square);
  free(data_positive);
  free(data_tmp);
  free(axis);
  free(d_axis);
}

int test_entity_basic_maths::run(){
  int err=TEST_PASSED;

  int res = integrator(data_square, size, d_axis);
  if(res!= data_square[size-1]) err |= TEST_WRONG_RESULT;
  res = integrator(data_positive, size, d_axis);
  if(res!= ((calc_type)(size*(size))/20.0 -1.0)) err |= TEST_WRONG_RESULT;
  //These slightly odd expression is because we're using basic trapezium numerical integral. We assume the top bnd can be projected on flat....

  memcpy((void*)data_square, (void*)data_tmp, sizeof(calc_type)*size);

  inplace_boxcar_smooth(data_tmp, size, 2, 1);
  int total=0;
  for(int i=0;i<size; ++i){
    total += data_tmp[i];
  }
  std::cout<<total<<std::endl;
  if(total != 0.0) err |=TEST_WRONG_RESULT;

  test_bed->report_err(err);
  return err;

}


