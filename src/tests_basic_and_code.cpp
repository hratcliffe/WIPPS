//
//  tests_basic_and_code.cpp
//  
//
//  Breakout by Heather Ratcliffe on 27/02/2017.
//
//

#ifdef RUN_TESTS_AND_EXIT

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <functional>
#include "tests.h"
#include "tests_basic_and_code.h"
#include "reader.h"
#include "data_array.h"

#include <math.h>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more

test_entity_reader::test_entity_reader(){
/** \brief Setup reader test*/
  name = "reader class";
  char block_id[ID_SIZE]= "run_info";
  test_rdr = new reader("./files/test", block_id);
  char block_id2[ID_SIZE] = "ax";
  accum_reader = new reader("./files/accum", block_id2);
}
test_entity_reader::~test_entity_reader(){
/** \brief Teardown reader test*/
  delete test_rdr;
  delete accum_reader;
}

int test_entity_reader::run(){
/**\brief Reads test data and checks
*
*For normal data we just check a test file opens and find the size from the final SDF block and test against on-disk size.  Then we read a distrib block and check it worked etc. For accumulated data we use test data with values corresponding to row in file, and check the values and correctness at each end.
@return Error code
\todo Investigate accumulator time being Inf in first bin
*/
  int err = TEST_PASSED;
  
  int sz = test_rdr->get_file_size();
  if(sz != size) err |= TEST_WRONG_RESULT;
  size_t n_dims;
  std::vector<size_t> dims;
  test_rdr->read_dims(n_dims, dims, "x_px/Background");
  //In test data this distrib is 2-D
  data_array distrib = data_array(dims[0], dims[1]);
  int tmp_err = test_rdr->read_distrib(distrib, "x_px/Background", 0);
  if(tmp_err || distrib.maxval() == distrib.minval()){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Error reading distribs", 1);
  }
  dims.clear();
  int rd_err = accum_reader->read_dims(n_dims, dims);
  if(rd_err){
    //Nothing else to do now
    err |=TEST_NULL_RESULT;
    test_bed->report_info("Error reading test files", 1);
    return err;
  }
  rd_err = 0;

  //Read some accumulated data and check it matches known values
  size_t time[3] = {0,10, 22};
  size_t space[2] = {0, dims[0]};
  data_array dat = data_array(dims[0], 40);

  rd_err = accum_reader->read_data(dat, time, space);
  test_bed->report_info("Reader returned "+mk_str(rd_err), 2);
  //rd_err of 2 means early stop but that is not a problem here
  if(!accum_reader->has_accum_data()){
    test_bed->report_info("Accumulated data check incorrect", 2);
    err |= TEST_ASSERT_FAIL;
  }
  if(rd_err!=1){
    //Test data should have row0=0, row1=2 (generation procedure does this, row3=1, 4=2 through to 9, then begin again at 1 through 10, 1 through 10.
    int row_vals[22] = {0, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1};
    int tot_errs = 0;
    for(size_t j = 0; j< dat.get_dims(1); j++){
      for(size_t i = 0; i< dat.get_dims(0); i++){
        tot_errs += (dat.get_element(i, j) != row_vals[j]);
      }
    }
    //Check axis values, first space
    for(size_t i = 0; i < dat.get_dims(0); i++){
      tot_errs += (dat.get_axis_element(0, i) != 1250*i);
    }
    //Now time, using known values
    float time_vals[13] = {0, 1.98054e-06, 0.000154482 ,0.000305003, 0.000455523, 0.000606044, 0.000756565, 0.000907086, 0.00105761, 0.00120813, 0.00135865,0.00150917, 0.00165969};
    for(size_t i = 1; i < 12; i++){
      if(std::abs(dat.get_axis_element(1, i) - time_vals[i]) > GEN_PRECISION) tot_errs++;
    }
    
    if(tot_errs != 0){
      err |=TEST_WRONG_RESULT;
      test_bed->report_info("Error reading accumulated data", 1);
    }
  }else{
    err |=TEST_NULL_RESULT;
    test_bed->report_info("Error reading test files", 1);
  }
  return err;
}

//----------------------------------------------------------------
test_entity_data_array::test_entity_data_array(){
/** \brief Setup data array test*/
  name = "data array class";
}

bool compare_2d(data_array const &lhs, data_array const &rhs, bool no_dims_match){
/** \brief Compare arrays
*
*Helper to compare two arrays in 2-d. Will early quit on first difference. NOTE the == op on arrays will only work for same dimensions
@param lhs First array
@param rhs Second array
@param no_dims_match If set, it will compare only the overlapping part for data
@return 0 for match, 1 for error
*/
  if(!no_dims_match){
    if(lhs.get_dims() != rhs.get_dims()) return 1;
    for(size_t i=0; i< lhs.get_dims(); i++) if(lhs.get_dims(i) != rhs.get_dims(i)) return 1;
  }
  for(size_t i=0; i<std::min(lhs.get_dims(0), rhs.get_dims(0)); i++){
    for(size_t j =0; j<std::min(lhs.get_dims(1), rhs.get_dims(1)); j++){
      if(lhs.get_element(i,j) != rhs.get_element(i, j)) return 1;
    }
  }
  for(size_t i=0; i<std::min(lhs.get_dims(), rhs.get_dims()); i++){
    for(size_t j=0; j< std::min(lhs.get_dims(i), rhs.get_dims(i)); j++){
      if(lhs.get_axis_element(i,j) != rhs.get_axis_element(i, j)) return 1;
    }
  }
  if(lhs.check_ids(rhs)) return 1;
  return 0;

}
bool compare_3d(data_array &lhs, data_array &rhs, bool no_dims_match){
/** \brief Compare arrays
*
*Helper to compare two arrays in 2-d. Will early quit on first difference. NOTE the == op on arrays will only work for same dimensions
@param lhs First array
@param rhs Second array
@param no_dims_match If set, it will compare only the overlapping part for data
@return 0 for match, 1 for error
*/
  if(!no_dims_match){
    if(lhs.get_dims() != rhs.get_dims()) return 1;
    for(size_t i=0; i< lhs.get_dims(); i++) if(lhs.get_dims(i) != rhs.get_dims(i)) return 1;
  }
  for(size_t i=0; i<std::min(lhs.get_dims(0), rhs.get_dims(0)); i++){
    for(size_t j =0; j<std::min(lhs.get_dims(1), rhs.get_dims(1)); j++){
      for(size_t k=0; k<std::min(lhs.get_dims(2), rhs.get_dims(2)); k++){
        if(lhs.get_element(i,j, k) != rhs.get_element(i, j, k)) return 1;
      }
    }
  }

  for(size_t i=0; i<std::min(lhs.get_dims(), rhs.get_dims()); i++){
    for(size_t j=0; j< std::min(lhs.get_dims(i), rhs.get_dims(i)); j++){
      if(lhs.get_axis_element(i,j) != rhs.get_axis_element(i, j)) return 1;
    }
  }
  if(lhs.check_ids(rhs)) return 1;

  return 0;

}

int test_entity_data_array::run(){
/** \brief Test data arrays
*
*Puts data into selected elements of a data array and reads back. Checks indexing, and bounds.
@return Error code
*/
  int err = TEST_PASSED;

  //Try tests in sections, aborting if we get a fatal error since this will only confuse the report
  err|=assign();
  if(test_bed->check_for_abort(err)) return err;
  err |= technical_tests();
  if(test_bed->check_for_abort(err)) return err;
  err |=basic_tests();
  if(test_bed->check_for_abort(err)) return err;
  err |= three_d_and_shift();
  if(test_bed->check_for_abort(err)) return err;
  err |= io_tests();
  return err;
}

int test_entity_data_array::set_vals_2d_and_sum(data_array &array, my_type &total){
/** \brief Create 2-d test array
*
*Assign each element of test_array and its axis elements to varying values. 
@param array The array to assign
@param[out] total Sum of values assigned
@return Error code
*/

  int err = TEST_PASSED;
  bool tmp_err;
  if(array.get_dims() != 2) return TEST_ASSERT_FAIL;
  total = 0.0; //Preserve total of all vals
  for(size_t i = 0; i < array.get_dims(0); i++){
    for(size_t j = 0; j < array.get_dims(1); j++){
      tmp_err = array.set_element(i, j, (i+1)*(2*j+1));
      if(tmp_err) err |= TEST_ASSERT_FAIL;
      total += (i+1)*(2*j+1);
    }
  }
  for(size_t i = 0; i < array.get_dims(); i++){
    for(size_t j = 0; j < array.get_dims(i); j++){
      tmp_err = array.set_axis_element(i, j, j*(i+1));
      if(tmp_err) err |= TEST_ASSERT_FAIL;
    }
  }
  return err;
}

int test_entity_data_array::set_vals_2d(data_array &array){
/** \brief Create 2-d test array
*
*Assign each element of test_array and its axis elements to varying values. 
@param array The array to assign. Should be 2-d of any size
@return Error code*/

  my_type tmp;
  return set_vals_2d_and_sum(array, tmp);
}

int test_entity_data_array::set_vals_3d(data_array &array){
/** \brief Create 3-d test array
*
*Assign each element of test_array and its axis elements to varying values.
@param array The array to assign. Should be 3-d of any size
@return Error code
*/

  int err = TEST_PASSED;
  bool tmp_err;
  if(array.get_dims() != 3) return TEST_ASSERT_FAIL;

  for(size_t i = 0; i<array.get_dims(0); i++){
    for(size_t j = 0; j<array.get_dims(1); j++){
      for(size_t k = 0; k<array.get_dims(2); k++){
        tmp_err=array.set_element(i, j, k, (i+1)*(2*j+1)*(4*k+1));
        if(tmp_err) err |= TEST_ASSERT_FAIL;
      }
    }
  }
  for(size_t i = 0; i<array.get_dims(); i++){
    for(size_t j = 0; j<array.get_dims(i); j++){
      tmp_err=array.set_axis_element(i, j, j*(i+1));
      if(tmp_err) err |= TEST_ASSERT_FAIL;
    }
  }
  return err;
}

int test_entity_data_array::assign(){
/** \brief Check basic assignments
*
*Set values and check basic assignment worked
@return Error code */

  int err = TEST_PASSED;
  my_type val;
  data_array test_array = data_array(10, 10);
  err |= this->set_vals_2d(test_array);

  //test assignments worked
  for(size_t i = 0; i<test_array.get_dims(0); i++){
    for(size_t j = 0; j<test_array.get_dims(1); j++){
      val = test_array.get_element(i,j);
      if(val != (i+1)*(2*j+1)) err |=TEST_WRONG_RESULT;
    }
  }
  for(size_t i = 0; i<test_array.get_dims(); i++){
    for(size_t j = 0; j<test_array.get_dims(i); j++){
      val=test_array.get_axis_element(i, j);
      if(val != (i+1)*j) err |=TEST_WRONG_RESULT;
    }
  }
  if(err == TEST_PASSED) test_bed->report_info("Assignment OK", 2);
  return err;
}

int test_entity_data_array::basic_tests(){
/** \brief test basic array functions
*
*Test maxval, total etc etc, resizer, maths
@return Error code */

  int err = TEST_PASSED;
  data_array test_array = data_array(10, 10);
  set_vals_2d(test_array);
  if(!test_array.is_good()) return TEST_ASSERT_FAIL;

  //test maxval function
  //Select a point withing the domain and bump it to max + 10 then check the new max is upped by 10 and is in the right place
  size_t i0 = test_array.get_dims(0)/2, i1 = test_array.get_dims(1)/3;
  my_type current_max = test_array.maxval();
  test_array.set_element(i0, i1, current_max + 10);
  std::vector<size_t> pos;
  my_type new_max = test_array.maxval(pos);
  if(new_max != current_max + 10 || pos.size()!=2||pos[0]!=i0||pos[1]!=i1){
    err |= TEST_WRONG_RESULT;
  }

  //test resizer
  size_t new2 = 6, new1 = 7;
  data_array old_array = test_array;

  test_bed->set_colour('*');
  test_bed->report_info("Testing resizer. Suggest using valgrind for memory checks", 2);
  test_bed->set_colour(0);

  test_bed->report_info("Initial size is "+mk_str(test_array.get_dims(0))+" by "+mk_str(test_array.get_dims(1)), 1);

  //Resize and check new sizes are as expected
  test_array.resize(1, new2, true);
  if(test_array.get_dims(1) !=new2){
    err |=TEST_WRONG_RESULT;
    test_bed->report_info("Second dim size is "+mk_str(test_array.get_dims(1))+" not "+mk_str(new2), 1);
  }
  test_array.resize(0, new1, true);
  if(test_array.get_dims(0) !=new1){
    err |=TEST_WRONG_RESULT;
    test_bed->report_info("First dim size is "+mk_str(test_array.get_dims(0))+" not "+mk_str(new1), 1);
  }
  //Compare old and new and report
  if(compare_2d(old_array, test_array, true)){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Resizer error, wrong values read", 1);
  }
  
  //Test total and averagers
  //Resize up and refill
  test_array = data_array(10, 10);
  my_type total_all;
  set_vals_2d_and_sum(test_array, total_all);
  //Calculate assigned total
  //Total up in both dimensions and compare to stored total
  data_array test_array2 = test_array.total(1);
  test_array2 = test_array2.total(0);
  my_type tot = test_array2.get_element((size_t)0);
  if(total_all != tot){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Totaler error", 1);
  }else{
    test_bed->report_info("Totaler OK", 2);
  }
  
  test_array = data_array(10, 10);
  set_vals_2d_and_sum(test_array, total_all);
  test_array2 = test_array.average(1);
  test_array2 = test_array2.average(0);
  tot = test_array2.get_element((size_t)0);
  if(total_all/100.0 != tot){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Averager error", 1);
  }else{
    test_bed->report_info("Averager OK", 2);
  }
  
  //resize to 1-d and test averaging fn
  test_array.resize(1, 1, true);
  my_type av = avval(test_array);
  my_type expected_av = ((10-1)/2.0 + 1);
  //Average of (i+1) for i=0 to 10
  if(av != expected_av){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Averager error in avval", 1);
  }else{
    test_bed->report_info("Averager OK", 1);
  }

  test_array = data_array(10, 10);
  set_vals_2d(test_array);

  //Test element-wise apply with simple +1 and log of constant values
  std::function<calc_type(calc_type)> plus1_function = [](calc_type el) -> calc_type { return el+1.0; } ;

  av = avval(test_array);
  test_array.apply(plus1_function);
  expected_av = av + 1.0;
  av = avval(test_array);
  if(av != expected_av){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Apply function error", 1);
  }
  std::function<calc_type(calc_type)> log_function = [](calc_type el) -> calc_type { return log(el); } ;
  
  data_array test_array_hand_logged = test_array;
  for(size_t i=0; i< test_array.get_dims(0); i++){
    for(size_t j=0; j< test_array.get_dims(1); j++){
      test_array_hand_logged.set_element(i, j, log(test_array_hand_logged.get_element(i, j)));
    }
  }
  test_array.apply(log_function);
  if(test_array != test_array_hand_logged){
  
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Apply function error", 1);
  }
  //Test apply from one array to another
  test_array2.clone_empty(test_array);
  expected_av = avval(test_array) + 1.0;
  test_array2.apply(plus1_function, test_array);
  av = avval(test_array2);
  if(av != expected_av){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Apply with array error", 1);
  }
  
  test_array = data_array(10, 10);
  set_vals_2d(test_array);
  
  //Test apply cross array using subtracting
  //Difference identical arrays and compare to 0-array. Also tests zero_data function
  test_array2 = test_array;
  data_array empty_array = test_array2;
  empty_array.zero_data();
  test_array2.apply(subtract, test_array);
  
  if(test_array2 != empty_array){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Subtractor error on first test", 1);
  }
  //Inverse check, this should fail to match, else our subtraction is nullifying stuff
  test_array2 = test_array;
  test_array2.set_element(test_array2.get_dims(0)/2, (size_t)0, -23.0);
  test_array2.apply(subtract, test_array);
  
  if(test_array2 == empty_array){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Subtractor error on second test", 1);
  }
  
  if(err == TEST_PASSED) test_bed->report_info("Basic tests OK");
  return err;
}

int test_entity_data_array::three_d_and_shift(){
/** \brief Test 3-d resizer and shift
*
* Tests resize in 3-d and shift in 3-d. If latter works, it should work in 2 or 1-d definitely.
@return Error code */

  int err = TEST_PASSED;
  //And now a 3-d version
  data_array test_array = data_array(10, 10, 10);
  err |= set_vals_3d(test_array);
  size_t new3 = 5;
  size_t new2 = 6;

  test_bed->report_info("Checking 3d",1);
  data_array old_array = test_array;

  test_array.resize(1, new2, true);
  test_array.resize(2, new3, true);

  if(test_array.get_dims(2) !=new3){
    err |=TEST_WRONG_RESULT;
    test_bed->report_info("Third dim size is "+mk_str(test_array.get_dims(2))+" not "+mk_str(new3), 1);
  }
  if(test_array.get_dims(1) !=new2){
    err |=TEST_WRONG_RESULT;
    test_bed->report_info("Second dim size is "+mk_str(test_array.get_dims(1))+" not "+mk_str(new2), 1);
  }

  if(compare_3d(old_array, test_array, true)){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Resizer error, wrong values read", 1);
  }

  //Now test the shift function by shifting, shifting back and then comparing
  //Try on all 3 dims
  old_array = test_array;
  test_array.shift(1, 3);
  test_array.shift(1, -3);

  if(old_array != test_array){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }

  test_array.shift(2, 3);
  test_array.shift(2, -3);

  if(old_array != test_array){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }

  test_array.shift(0, 2);

  if(old_array == test_array){
  //If they still compare equal we has problem
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, no shift applied", 1);
  }

  test_array.shift(0, -2);
  
  if(old_array != test_array){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }
  return err;

}

int test_entity_data_array::technical_tests(){
/** \brief Technical testing
*
*Check things like copy constructors
@return Error code */

  test_bed->report_info("Checking technical aspects", 2);

  //Create empty array
  data_array test_array2 = data_array(10, 10);
  set_vals_2d(test_array2);

  int err = TEST_PASSED;
  data_array dat = test_array2;
  if(test_array2 != dat){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Copy or equality problem", 1);
  }
  try{
    std::vector<data_array> my_vec;
    my_vec.push_back(test_array2);
    my_vec.push_back(test_array2);
    my_vec.resize(100);
    my_vec.resize(10);
    if(my_vec[0] != test_array2) err|= TEST_WRONG_RESULT;
  }catch(const std::exception& e){
    //Swallow and continue if possible,
    std::string message = e.what();
    test_bed->report_info("Exception message " +message, 1);
    err |= TEST_ASSERT_FAIL;
  }
  
  //Check the downcasting conversions
  my_array source_array = my_array(10, 10);
  source_array.set_element((size_t)0, (size_t)0, 3.0);
  dat = source_array;
  my_type tmp = avval(source_array), tmp2 = avval(dat);

  if(dat != source_array || tmp != tmp2){
    test_bed->report_info("Error in conversion", 2);
    err |= TEST_ASSERT_FAIL;
  }else{
    test_bed->report_info("Conversion OK", 2);
  }

  if(err == TEST_PASSED) test_bed->report_info("Technical aspects OK", 1);
  return err;
}

int test_entity_data_array::io_tests(){
/** \brief Check read/write
*
*Check data_array write to file and then read
@return Error code */

  test_bed->report_info("Checking file io", 1);
  int err = TEST_PASSED;
  bool err2 = false;

  data_array test_array = data_array(10, 10);
  set_vals_2d(test_array);

  std::string filename = tests_tmp_dir+"test_file.dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    err2 = test_array.write_to_file(file);
    file.close();
  }else{
    err2 = true;
  }
  if(err2) test_bed->report_info("Error writing testfile", 1);
  
  //Setup the size etc
  data_array new_array = test_array;
  //Change some element to prove reading is doing something
  new_array.set_element(5, 5, -20.0);
  file.open(filename.c_str(),std::ios::in|std::ios::binary);
  err2 = new_array.read_from_file(file);
  if(err2) test_bed->report_info("Error reading testfile", 1);
  file.close();
  
  if(err2) err|=TEST_ASSERT_FAIL;
  if(test_array != new_array) err |= TEST_WRONG_RESULT;

  return err;
}

//----------------------------------------------------------------

test_entity_get_and_fft::test_entity_get_and_fft(){
/** \brief Setup FFT testing */
  name = "read and FFT";
}

test_entity_get_and_fft::~test_entity_get_and_fft(){
/** \brief Teardown FFT testing */
  if(test_rdr) delete test_rdr;
}

int test_entity_get_and_fft::run(){
/** \brief Checks full read and fft procedure
*
*Reads a test sdf file, stores into data array and runs fft. Test data should be a sine curve with one major frequency which is then checked. Note frequency is hard coded to match that produced by ./files/sin.deck
@return Error code
*/

  int err = TEST_PASSED;
  
  std::string block = "ex";
  test_rdr = new reader("./files/sin", block.c_str());
  err |= one_d();

  if(test_rdr) delete test_rdr;
  block = "ax";
  test_rdr = new reader("./files/sinAcc", block.c_str());
  err |= two_d();

  return err;
}

int test_entity_get_and_fft::one_d(){
/** \brief Checks fft in 1-d
*
*Read, fft and check for 1-d data
@return Error code
*/

  int err = TEST_PASSED;

  size_t tim_in[3] = {0, 1, 0};
  size_t space_in[2] = {0, 0};
  int n_tims = std::max((int) (tim_in[1]-tim_in[0]), 1);

  //Get the dimensions from file
  size_t n_dims;
  std::vector<size_t> dims;
  bool rd_err = test_rdr->read_dims(n_dims, dims);

  if(rd_err){
    err |= TEST_ASSERT_FAIL;
    test_bed->report_info("Error reading file or block", 1);
    return err;
  }
  space_in[1] = dims[0];

  if(n_dims != 1){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Wrong array dims read", 1);
    return err;
  }

  //Assign data arrays
  data_array test_dat = data_array(dims[0], n_tims);
  data_array test_dat_fft = data_array(dims[0], n_tims);
  if(!test_dat.is_good()||!test_dat_fft.is_good()){
    err |= TEST_ASSERT_FAIL;
    return err;
  }
  
  //Read the data and check frequencies
  test_rdr->read_data(test_dat, tim_in, space_in);

  my_type expected_max = 1.2566371e-4;//Expected frequency of max from test data

  //Check with raw data. This we know is exactly symmetric
  err |= fft_and_check_1d(test_dat, test_dat_fft, expected_max);
  
  //Now size down by 1 element and redo. This checks odd and even total sizes in the FFT
  test_dat.resize(0, dims[0]-1, true);
  test_dat_fft.resize(0, dims[0]-1, true);
  err |= fft_and_check_1d(test_dat, test_dat_fft, expected_max);
  
  return err;

}

int test_entity_get_and_fft::fft_and_check_1d(data_array & dat_in, data_array & dat_fft, my_type expected_max, bool single_max){
/**  \brief Check 1-d FFT frequencies
*
* FFTS data_in into data_fft and then hunts for 1 or two maxima, depending on n_maxima, and checks they're at axis values of ± expected_max
@param dat_in Input array 
@param dat_fft FFT'd data array
@param expected_max Expected peak axis value
@param single_max Whether data has a single maximum freq (e.g. at 0)
@return Error code
*/
  int err = TEST_PASSED;
  bool tmp_err = fft_array(dat_in, dat_fft);
  if(tmp_err) err |= TEST_ASSERT_FAIL;
  if(dat_fft.check_ids(dat_in)) err |= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("1D FFT performed without error", 1);

  //Get the two maxes (-ve and +ve frequency) and check position is correct
  int max_index = 0;
  my_type max_val = 0;
  std::vector<size_t> max_pos;
  bool both_freqs_correct = true;
  int sgn = 1;

  //Cheating repeat of checker code either once if single_max or twice otherwise for +ve and -ve frequencies
  for(int i=0; i< 2-single_max; i++){
    if(i==0) max_val = dat_fft.maxval(max_pos);//NB this is updating max_pos
    //Find max above previous max_index
    else max_val = dat_fft.maxval(max_pos, max_index+1);

    if(max_pos.size() < 1){
      err |= TEST_WRONG_RESULT;
      continue;
    }
    //Since we know we're 1-d we consider only 0 axis
    max_index = max_pos[0];
    sgn = dat_fft.get_axis_element(0,max_index)/std::abs(dat_fft.get_axis_element(0,max_index));
    //Allow either exact match or boxed match in general. First implies second
    //Exact match to within PRECISION
    bool exact_match = !(std::abs(std::abs(dat_fft.get_axis_element(0, max_index)) - expected_max) > PRECISION);
    //Boxed match, found axis elements box the expected value
    bool boxed_match = ((dat_fft.get_axis_element(0, max_index+1) >= sgn*expected_max) && (dat_fft.get_axis_element(0, max_index) <= sgn*expected_max)) || ((dat_fft.get_axis_element(0, max_index) >= sgn*expected_max) && (dat_fft.get_axis_element(0, max_index - 1) <= sgn*expected_max));
    if(!exact_match && !boxed_match){
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("Max freq is "+mk_str(dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+")", 1);
      both_freqs_correct = false;
    }
  }
  
  if(both_freqs_correct){
    test_bed->report_info("1D FFT Frequency correct!", 1);
  }else{
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Error wrong 1D FFT Frequency", 1);
  }

  return err;
}

int test_entity_get_and_fft::two_d(){
/**  \brief Check 2-d FFT frequencies
*
* FFTS data_in into data_fft and then hunts for maxima, and checks against hard coded expected values
@return Error code
*/

  int err = TEST_PASSED;

  size_t tim_in[3] = {0, 3, 100};
  size_t space_in[2] = {0, 0};

  int n_tims = tim_in[2];

  size_t n_dims;
  std::vector<size_t> dims;
  bool rd_err = test_rdr->read_dims(n_dims, dims);

  if(rd_err){
    err |= TEST_ASSERT_FAIL;
    test_bed->report_info("Error reading file or block", 1);
    return err;
  }
  space_in[1] = dims[0];
  int space_size = space_in[1]-space_in[0];
  //This is accumulated data and time dim is omitted by get_dims
  if(n_dims !=1){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Array dims wrong", 1);
    return err;
    //nothing more worth doing right now...
  }

  data_array test_dat = data_array(space_size, n_tims);
  data_array test_dat_fft = data_array(space_size, n_tims);
  if(!test_dat.is_good()||!test_dat_fft.is_good()){
    err|=TEST_ASSERT_FAIL;
    return err;
  }
  
  test_rdr->read_data(test_dat, tim_in, space_in);

  bool tmp_err = fft_array(test_dat,test_dat_fft);
  if(tmp_err) err|=TEST_ASSERT_FAIL;
  if(test_dat_fft.check_ids(test_dat)) err |= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("2D read and FFT reports no error", 1);
  
  int max_index = 0;
  my_type max_val = 0;
  std::vector<size_t> max_pos;
  my_type expected_max = 1.2566371e-4;
  bool both_freqs_correct = true;

  max_val = test_dat_fft.maxval(max_pos);

  if(max_pos.size() < 2){
    err |= TEST_WRONG_RESULT;
    max_pos.resize(2);//Force resize so we get the normal reports etc
  }
  max_index = max_pos[0];
  
  if(std::abs(std::abs(test_dat_fft.get_axis_element(0,max_index)) - expected_max) > PRECISION){
    err|= TEST_WRONG_RESULT;
    test_bed->report_info("Max wavenum is "+mk_str(test_dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+", "+mk_str(max_pos[1])+")", 1);
    both_freqs_correct = false;
  }

  if(max_pos[1] != test_dat_fft.get_dims(1)/2){
    err|= TEST_WRONG_RESULT;
    test_bed->report_info("Max freq is "+mk_str(test_dat_fft.get_axis_element(1,max_pos[1]))+" ("+mk_str(max_index)+", "+mk_str(max_pos[1])+")", 1);
    both_freqs_correct = false;
  }

  std::vector<std::pair<size_t, size_t> > ranges;
  ranges.push_back(std::make_pair<size_t, size_t>(test_dat_fft.get_dims(0)/2, test_dat_fft.get_dims(0)));
  ranges.push_back(std::make_pair<size_t, size_t>(0, test_dat_fft.get_dims(1)));
  //Now check in upper half. I know the test data is assymetric. 

  max_val = test_dat_fft.partial_maxval(ranges, max_pos);
  if(std::abs(std::abs(test_dat_fft.get_axis_element(0,max_pos[0])) - expected_max) > PRECISION){
    err|= TEST_WRONG_RESULT;
    test_bed->report_info("Max wavenum is "+mk_str(test_dat_fft.get_axis_element(0,max_pos[0]))+" ("+mk_str(max_pos[0])+", "+mk_str(max_pos[1])+")", 1);
    both_freqs_correct = false;
  }

  
  if(both_freqs_correct) test_bed->report_info("FFT Frequency correct!", 1);
  
  std::string filename, time_str;
  int err2;
  time_str = mk_str(test_dat_fft.time[0], true)+"_"+mk_str(test_dat_fft.time[1],true);
  std::string block = test_dat_fft.block_id;
  
  filename = tests_tmp_dir+"FFT2DTest_"+block +"_"+time_str+"_"+mk_str(test_dat_fft.space[0])+"_"+mk_str(test_dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
      err2 = test_dat_fft.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();

  return err;

}
//----------------------------------------------------------------

test_entity_basic_maths::test_entity_basic_maths(){
/** \brief Setup basic maths testing */
  name = "basic maths helpers";
}
test_entity_basic_maths::~test_entity_basic_maths(){
/** \brief Teardown basic maths tests*/
  teardown_arrays();
}

void test_entity_basic_maths::setup_arrays(){
/** \brief Allocate and fill test arrays
*
*Sets up test arrays with special data for testing
*/

  data_square=(calc_type*)calloc(size,sizeof(calc_type));
  data_positive=(calc_type*)calloc(size,sizeof(calc_type));
  data_tmp=(calc_type*)calloc(size,sizeof(calc_type));
  axis=(calc_type*)calloc(size,sizeof(calc_type));
  d_axis=(calc_type*)calloc(size,sizeof(calc_type));

  axisf=(my_type*)calloc(size,sizeof(my_type));

  data_square[0] = 1.0;
  data_positive[0] = 0.0;
  axis[0] = 0.0;
  axisf[0] = 0.0;
  d_axis[0] = 1.0;
  for(int i=1; i<size; ++i){
    data_square[i] = - data_square[i-1];
    //alternating square wave, average=int=0
    data_positive[i] = data_positive[i-1] + 0.1;
    //monotonic increase, integral = size * size/10/2 depending on upper bnd
    d_axis[i] = 1.0;
    axis[i] = axis[i-1] + d_axis[i];
    axisf[i] = axisf[i-1] + 1.0;

  }

}
void test_entity_basic_maths::teardown_arrays(){
/** \brief Free test arrays
*
*Free arrays setup by test_entity_basic_maths::setup_arrays()
*/
  free(data_square);
  free(data_positive);
  free(data_tmp);
  free(axis);
  free(d_axis);
  free(axisf);
}

int test_entity_basic_maths::run(){
/** \brief Test basic maths functions
*
*Tests functions such as "where" value locator, integrator, smoothing, interpolator, cubic equation solver and data slice flattener
@return Error code
*/

  int err = TEST_PASSED;
  setup_arrays();
  my_type target;
  int whe;
  //Test where location
  target = 13.5;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = 254.89;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = -2;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = 1302;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  if(err == TEST_PASSED) test_bed->report_info("Where OK", 1);

  //Test integration
  calc_type res = integrator(data_square, size, d_axis);
  if(res!= 0.0) err |= TEST_WRONG_RESULT;
  res = integrator(data_positive, size, d_axis);
  if(std::abs(res - (pow((calc_type)(size-1), 2)/20.0)) > res*PRECISION) err |= TEST_WRONG_RESULT;
  //test it's correct to within some finite precision, defined in header
  if(err == TEST_PASSED) test_bed->report_info("Integrator OK", 1);

  //Test in-place smoothing
  memcpy((void*)data_tmp, (void*)data_square, sizeof(calc_type)*size);

  inplace_boxcar_smooth(data_tmp, size, 2, 1);
  calc_type total=0;
  for(int i=1;i<size-1; ++i){
    total += std::abs(data_tmp[i]);
  }
  inplace_boxcar_smooth(data_tmp+1, size-2, 2, 1);
  for(int i=2;i<size-3; ++i){
    total += std::abs(data_tmp[i]);
  }
  if(std::abs(total) > PRECISION) err |=TEST_WRONG_RESULT;
  //Smooth of 2 on square wave should give 0. Do with and without offset as extra check. Omit ends in total
  
  memcpy((void*)data_tmp, (void*)data_positive, sizeof(calc_type)*size);
  
  inplace_boxcar_smooth(data_tmp, size, 3, 0);
  total=0;
  for(int i=4;i<size-4; ++i){
    total += std::abs(data_positive[i] - data_tmp[i]);
  }
  //Smooth on straight line with odd width should do nothing except at ends...
  if(std::abs(total) > PRECISION) err |=TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Boxcar smooth OK", 1);

  //Test interpolation on tiny arrays with known values
  {
    calc_type axis[4]={0.0,1.0,2.0,3.0}, vals[4]={0.0,1.0,0.0,1.0}, target, interp;

    target = 0.25;
    interp = interpolate_nearest(axis, vals, target);
    if(std::abs(interp - 0.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 0.75;
    interp = interpolate_nearest(axis, vals, target);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    //Nearest value interpol
    target = 0.5;
    interp = interpolate_linear(axis, vals, target);
    if(std::abs(interp - 0.5)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.0;
    interp = interpolate_linear(axis, vals, target);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.5;
    interp = interpolate_linear(axis+1, vals+1, target);
    if(std::abs(interp - 0.5)> PRECISION) err |=TEST_WRONG_RESULT;
    //2 pt linear interpol
  }
  
  if(err == TEST_PASSED) test_bed->report_info("Interpolate OK", 1);

  //  x^3 - 17x^2 + 92x - 150.
  std::vector<calc_type> result = cubic_solve(-17.0, 92.0, -150.0);
  if((result.size() != 1) || std::abs(result[0] - 3.0) > PRECISION) err|=TEST_WRONG_RESULT;
  //Test with polynomial of known integer coefficients and roots

  result = cubic_solve(-20.5, 100.0, -112.76);
  //test with random polynomial which happens to have 3 real roots
  calc_type res_el, tot;
  for(size_t i=0; i<result.size(); ++i){
    res_el = result[i];
    tot = std::pow(res_el, 3) -20.5*std::pow(res_el, 2) +100.0*res_el - 112.76;
    if(std::abs(tot) > PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Cubic root does not solve polynomial, mismatch "+mk_str(tot)+" for root "+mk_str(res_el), 2);
      
      test_bed->report_info("Vieta a) "+ mk_str(result[0]+result[1]+result[2] - 20.5), 2);
      test_bed->report_info("Vieta b) "+mk_str( result[0]*result[1]+result[0]*result[2]+result[2]*result[1] - 100.0), 2);
      test_bed->report_info("Vieta c) "+mk_str(result[0]*result[1]*result[2]-112.76), 2);

    }
  }
 	
  /*"I think you should be able to recognize them using Vieta's formula for cubic equations, which states that if a cubic equation x3+ax2+bx+c=0x3+ax2+bx+c=0 has three different roots x1,x2,x3x1,x2,x3, then:
  −a=x1+x2+x3 b = x1x2+x1x3+x2x3 and -c = x1x2x3"
**/

  //Check our array slice flattener
{  size_t dims[3];
  size_t n_dims = 3;
  dims[0] = 5; dims[1] = 6; dims[2] = 5;
  size_t total_sz = dims[0]*dims[1]*dims[2];
  my_type * in, *out;
  in = (my_type *) calloc(total_sz, sizeof(my_type));
  out = (my_type *) calloc(total_sz/dims[1], sizeof(my_type));
  
  int tot_on_dim1 = 0;
  for(size_t j=0; j<dims[1]; j++) tot_on_dim1 += j;
  
  for(size_t i = 0; i< dims[0]; i++){
    for(size_t j=0; j<dims[1]; j++){
      for(size_t k=0; k< dims[2]; k++){
        *(in+(k*dims[1]+ j)*dims[0] + i) = i+j;
        
      }
    }
  }
  flatten_fortran_slice(in, out, n_dims, dims, 1);
  int errs = 0;
  for(size_t i = 0; i< dims[0]; i++){
    for(size_t k=0; k< dims[2]; k++){
      if(*(out + (k*dims[0] + i)) != tot_on_dim1+ dims[1]*i) errs++;
    }
  }
  if(errs > 0) err |= TEST_WRONG_RESULT;

  free(in);
  free(out);
  
}

{  size_t dims[2];
  size_t n_dims = 2;
  dims[0] = 5; dims[1] = 6;
  size_t total_sz = dims[0]*dims[1];
  my_type * in, *out;
  in = (my_type *) calloc(total_sz, sizeof(my_type));
  out = (my_type *) calloc(total_sz/dims[1], sizeof(my_type));
  
  int tot_on_dim1 = 0;
  for(size_t j=0; j<dims[1]; j++) tot_on_dim1 += j;
  
  for(size_t i = 0; i< dims[0]; i++){
    for(size_t j=0; j<dims[1]; j++){
        *(in+ j*dims[0] + i) = i+j;
    }
  }
  flatten_fortran_slice(in, out, n_dims, dims, 1);
  int errs = 0;
  for(size_t i = 0; i< dims[0]; i++){
      if(*(out + i) != tot_on_dim1+ dims[1]*i) errs++;
  }
  
  if(errs > 0) err |= TEST_WRONG_RESULT;
  free(in);
  free(out);

}
  
  if(err == TEST_PASSED) test_bed->report_info("Flattener OK", 1);
  return err;

}

test_entity_extern_maths::test_entity_extern_maths(){
/** \brief Setup test for external maths*/
  name = "external maths routines";
}

int test_entity_extern_maths::run(){
/** \brief Test external maths
*
* Tests bessel functions, and short-cut approximations
@return Error code
*/
  //cyl_bessel_j(v, x) = Jv(x)
  //cyl_neumann(v, x) = Yv(x) = Nv(x)
  double bess, arg, bess1, bess2;
  int index, err=TEST_PASSED;

  index = 0;
  arg = 2.40482555769577;
  bess = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  index = 1;
  arg=7.01558666981561;
  bess = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  index=5;
  arg=12.3386041974669;
  bess = boost::math::cyl_bessel_j(index, arg);

  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  //Test identities used to save time
  index=3;
  arg=1.3457174;
  bess = boost::math::cyl_bessel_j(index-1, arg);
  bess1 = boost::math::cyl_bessel_j(index+1, arg);
  bess2 = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;

  index=7;
  arg=-20.98;
  bess = boost::math::cyl_bessel_j(index-1, arg);
  bess1 = boost::math::cyl_bessel_j(index+1, arg);
  bess2 = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;


  index =0;
  arg =1.0;
  bess = boost::math::cyl_bessel_j(index, arg);

  if(std::abs(bess - 0.7651976865579665514497) >PRECISION) err|= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Bessel functions OK", 1);

  return err;
}

#endif
