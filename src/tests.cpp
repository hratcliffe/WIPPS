//
//  tests.cpp
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

#ifdef RUN_TESTS_AND_EXIT

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include "tests.h"
#include "reader.h"
#include "main_support.h"
#include "plasma.h"
#include "controller.h"
#include "data_array.h"
#include "spectrum.h"
#include "d_coeff.h"
#include "non_thermal.h"

#include <math.h>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
/** \todo Deliberately failing tests to check not doing nothing \todo Add scripts directory for level1 \todo Add tests directory for easy cleanup*/

extern tests * test_bed; /**< Global testbed, define somewhere in your code*/

extern const mpi_info_struc mpi_info;
extern deck_constants my_const;

const int err_codes[err_tot] ={TEST_PASSED, TEST_WRONG_RESULT, TEST_NULL_RESULT, TEST_ASSERT_FAIL, TEST_USERDEF_ERR1, TEST_USERDEF_ERR2, TEST_USERDEF_ERR3, TEST_USERDEF_ERR4, TEST_FATAL_ERR};/**< List of error codes available*/

std::string err_names[err_tot]={"None", "Wrong result", "Invalid Null result", "Assignment or assertion failed", "{Message 1 here}", "{Message 2 here}", "{Message 3 here}", "{Message 4 here}", "Fatal error"};/**< Names corresponding to error codes, which are reported in log files*/

tests::tests(){
  setup_tests();
}
tests::~tests(){
  cleanup_tests();
}

void tests::set_verbosity(int verb){
/** Set the verbosity of testing output, from 0 (minimal) to max_verbos.*/
  if((verb > 0)) this->verbosity = std::max(verb, max_verbos);

}

void tests::setup_tests(){
/** \brief Setup test bed
*
*Opens reporting file. Then instantiates all the test objects and adds them into the test_list
*/
  outfile = new std::fstream();
  outfile->open(filename.c_str(), std::ios::out);
  if(!outfile->is_open()){
    my_print("Error opening "+filename, mpi_info.rank);
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
  test_obj = new test_entity_extern_maths();
  add_test(test_obj);
  test_obj = new test_entity_plasma();
  add_test(test_obj);
  test_obj = new test_entity_spectrum();
  add_test(test_obj);
/*  test_obj = new test_entity_levelone();
  add_test(test_obj);
  test_obj = new test_entity_d();
  add_test(test_obj);*/
  test_obj = new test_entity_nonthermal();
  add_test(test_obj);

}

void tests::add_test(test_entity * test){
  /** Adds a test to the list to be performed*/
  test_list.push_back(test);
}

bool tests::is_fatal(int err){

  if((err & TEST_FATAL_ERR) == TEST_FATAL_ERR) return 1;
  else return 0;
}

bool tests::check_for_abort(int err){

  if(is_fatal(err)){
    set_colour('r');
    set_colour('*');
    my_print("Fatal error occured. Aborting test "+test_list[current_test_id]->name, mpi_info.rank);
    set_colour();
    return true;
  }
  else{
    return false;
  }
}
void tests::report_err(int err, int test_id){
/** \brief Log error
*
* Logs error text corresponding to code err for test defined by test_id. Errors are always recorded.*/
  if(test_id == -1) test_id = current_test_id;
  if(err ==TEST_PASSED) set_colour('b');
  else set_colour('r');
  if(is_fatal(err)) set_colour('*');
  my_print(outfile, get_printable_error(err, test_id), mpi_info.rank);
  my_print(nullptr, get_printable_error(err, test_id), mpi_info.rank);
  set_colour();

}

void tests::report_info(std::string info, int verb_to_print, int test_id){
/** \brief Other test info
*
*Records string info to the tests.log file and to screen, according to requested verbosity.
*/
  if(test_id == -1) test_id = current_test_id;
  if(verb_to_print <= this->verbosity){
    my_print(outfile, info, mpi_info.rank);
    my_print(nullptr, info, mpi_info.rank);
  
  }

}

std::string tests::get_printable_error(int err, int test_id){
/** \brief Make an error message
*
* Converts error code to printable string, adds code for reference and adds test name. Note code is bitmask and additional errors are appended together
*/
  std::string err_string="";
  if(err!=TEST_PASSED){
    for(int i=err_tot-1; i>0; --i){
      //Run most to least significant
      if((err & err_codes[i]) == err_codes[i]){
        err_string +=err_names[i];
        if(i >1) err_string +=" &";
        err_string += " ";
      }
    }
  
    err_string = "Error "+err_string+"(code "+mk_str(err)+") on";
  }
  else err_string = "Passed";
  return err_string+" test "+test_list[test_id]->name;

}

void tests::set_colour(char col){
/** \brief Set output text colour
*
*Set terminal output colour using std escape sequences. Accepts no argument to return to default, or rgb, cmyk and white to test text colour. NB technically not MPI safe. Use sparingly to highlight important information.
*/

  my_print(this->get_color_escape(col), mpi_info.rank, 0, true);
}

inline std::string tests::get_color_escape(char col){
/** \brief
*\copydoc dummy_colour This returns the terminal escape string to set given colour.
*/
  if(col >='A' and col <='Z') col += 32;
  //ASCII upper to lower
  switch (col) {
    case 0:
    case '0':
      return "\033[0m";
      break;//Redundant but clearer
    case 'r':
      return "\033[31m";
      break;
    case 'g':
      return "\033[32m";
      break;
    case 'b':
      return "\033[34m";
      break;
    case 'c':
      return "\033[36m";
      break;
    case 'm':
      return "\033[35m";
      break;
    case 'y':
      return "\033[33m";
      break;
    case 'w':
      return "\033[37m";
      break;
    case 'k':
      return "\033[30m";
      break;
    case '*':
    //bold
      return "\033[1m";
      break;
    case '_':
    //underline
      return "\033[4m";
      break;
    case '?':
    //blink. very annoying
      return "\033[5m";
      break;
    case '$':
    //reverse fore/back ground
      return "\033[7m";
      break;
    
    default:
      return "";
  }

}

void tests::cleanup_tests(){
/** \brief Delete test objects
*
*
*/
  if(outfile->is_open()){
    this->report_info("Testing complete and logged in " +filename, 0);
    outfile->close();
  }else{
    this->report_info("No logfile generated", 0);

  }
  delete outfile;
//  test_obj * obj;
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
  //  obj = test_list[current_test_id];
   // test_list.erase(current_test_id);
    delete test_list[current_test_id];
    test_list[current_test_id] = nullptr;
  
  }

  
}

bool tests::run_tests(){
/** \brief Run scheduled tests
*
*Runs each test in list and reports total errors found
*/

  int total_errs = 0;
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
    total_errs += (bool) test_list[current_test_id]->run();
    //Add one if is any error returned
  }
  this->set_colour('r');
  this->report_info(mk_str(total_errs)+" failed tests", mpi_info.rank);
  this->set_colour();
  return total_errs > 0;

}
//----------------------------------------------------------------

test_entity_reader::test_entity_reader(){
  name = "reader class";
  char block_id[ID_SIZE]= "run_info";
  test_rdr = new reader("./files/test", block_id);
  char block_id2[ID_SIZE] = "ax";

  accum_reader = new reader("./files/accum", block_id2);
  
}
test_entity_reader::~test_entity_reader(){

  delete test_rdr;
  delete accum_reader;
}

int test_entity_reader::run(){
/**\brief Reads test data and checks
*
*For normal data we just check a test file opens and find the size from the final SDF block and test against on-disk size. For accumulated data we use test data with values corresponding to row in file, and check the values and correctness at each end.
*/
  int err = TEST_PASSED;
  
  int sz = test_rdr->get_file_size();
  if(sz != size) err |= TEST_WRONG_RESULT;
  int rd_err;
  size_t n_dims;
  std::vector<size_t> dims;
  rd_err = accum_reader->read_dims(n_dims, dims);
  int time[3]={0,10, 40}, space[2]={0, (int)dims[0]};

  if(!rd_err){
    
    data_array dat = data_array(dims[0], 40);
    //We know there are at most ten rows per accumulate in our test file and 4 files
    rd_err = accum_reader->read_data(dat, time, space);
    size_t a, b;
    if(rd_err!=1){
       //Test data should have row 0 =0, row 1=2 (generation procedure does this, row 3=1, 4=2 through to 9, then begin again at 1 through 10, 1 through 10.
       //Check some selection of elements and sum errs. Note all test entries are whole ints
       int tot_errs=0;

       //Check a few selected rows are right...
       a=0;
       b=0;
       tot_errs += !(dat.get_element(a,b) == 0.0);
       b=1;
       tot_errs += !(dat.get_element(a,b) == 2.0);
       b=2;
       tot_errs += !(dat.get_element(a,b) == 1.0);
       b=10;
       tot_errs += !(dat.get_element(a,b) == 9.0);
       b=11;
       tot_errs += !(dat.get_element(a,b) == 1.0);
      
       //Check both ends of a few rows for off by ones etc
       int rows[3] = {1,10,11};
       for(int i = 0; i<3; i++ ){
         tot_errs += !(dat.get_element(0,rows[i]) == dat.get_element(1,rows[i]));
         tot_errs += !(dat.get_element(0,rows[i]) == dat.get_element(399,rows[i]));
         tot_errs += !(dat.get_element(0,rows[i]) == dat.get_element(200,rows[i]));

       }
       
       if(tot_errs != 0){
         err |=TEST_WRONG_RESULT;
         test_bed->report_info("Error reading accumulated data", 1);
       }

    }else{
      err |=TEST_NULL_RESULT;
      test_bed->report_info("Error reading test files", 1);
    }
    
  }else{
    err |=TEST_NULL_RESULT;
    test_bed->report_info("Error reading test files", 1);

  }

  test_bed->report_err(err);
  return err;

}

//----------------------------------------------------------------
test_entity_data_array::test_entity_data_array(){
  name = "data array class";
  test_array = data_array(10, 10);
  
}
test_entity_data_array::~test_entity_data_array(){

//  delete test_array;
}

bool compare_2d(data_array &lhs, data_array &rhs, bool no_dims_match){
/**Helper to compare two arrays in 2-d. Will early quit on first difference. If no_dims_match is set, it will compare only the overlap @return 0 for match, 1 for error*/
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
/**Helper to compare two arrays in 3-d. Will early quit on first difference. If no_dims_match is set, it will compare only the overlap @return 0 for match, 1 for error*/
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
/** \brief Puts data into selected elements of a data array and reads back. Checks indexing, and bounds.
*
*
*/
  int err = TEST_PASSED;

  err|=assign();
  if(test_bed->check_for_abort(err)) return err;
  err |=basic_tests();
  if(test_bed->check_for_abort(err)) return err;
  err |= technical_tests();
  if(test_bed->check_for_abort(err)) return err;
  err |= three_d_and_shift();
  if(test_bed->check_for_abort(err)) return err;
  err |= io_tests();
  test_bed->report_err(err);
  return err;
}

int test_entity_data_array::assign(){
/** Set values and check basic assignment worked*/
  int err = TEST_PASSED;
  bool tmp_err;
  my_type val;
  //assign each element to unique val and axes to position

  for(size_t i=0; i<test_array.get_dims(0); i++){
    for(size_t j =0; j<test_array.get_dims(1); j++){
      tmp_err=test_array.set_element(i, j, (i+1)*(2*j+1));
      if(tmp_err) err |= TEST_ASSERT_FAIL;
    }
  }
  for(size_t i=0; i<test_array.get_dims(); i++){
    for(size_t j =0; j<test_array.get_dims(i); j++){
      tmp_err=test_array.set_axis_element(i, j, j*(i+1));
      if(tmp_err) err |= TEST_ASSERT_FAIL;
    }
  }

  //test assignments worked

  for(size_t i=0; i<test_array.get_dims(0); i++){
    for(size_t j =0; j<test_array.get_dims(1); j++){
      val = test_array.get_element(i,j);
      if(val != (i+1)*(2*j+1)) err |=TEST_WRONG_RESULT;
    }
  }
  for(size_t i=0; i<test_array.get_dims(); i++){
    for(size_t j =0; j<test_array.get_dims(i); j++){
      val=test_array.get_axis_element(i, j);
      if(val != (i+1)*j) err |=TEST_WRONG_RESULT;
    }
  }

  return err;
}

int test_entity_data_array::basic_tests(){
/** Test maxval, resizer, maths*/

  int err = TEST_PASSED;
  if(!test_array.is_good()) return TEST_ASSERT_FAIL;

  //test maxval function
  size_t i0=test_array.get_dims(0)/2, i1=test_array.get_dims(1)/3;
  my_type current_max = test_array.maxval();
  test_array.set_element(i0, i1, current_max+10);
  std::vector<size_t> pos;
  current_max -= test_array.maxval(pos);
  if(current_max != -10 || pos.size()!=2||pos[0]!=i0||pos[1]!=i1){
    err |= TEST_WRONG_RESULT;
  }

  //test resizer

  size_t new2=6, new1=7;
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
  
  //resize to 7x1 and test average fn
  test_array.resize(1, 1, true);
  my_type av = test_array.avval();
  my_type expected_av = ((new1-1)/2.0 + 1);
  //Average of (i+1) for i=0 to new1

  if(av != expected_av){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Averager error", 1);
  }
  
  //Test apply with simple +1 and log of constant values
  std::function<calc_type(calc_type)> plus1_function = [](calc_type el) -> calc_type { return el+1.0; } ;

  test_array.apply(plus1_function);
  expected_av = av + 1.0;
  av = test_array.avval();
  if(av != expected_av){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Apply function error", 1);
  }
  std::function<calc_type(calc_type)> log_function = [](calc_type el) -> calc_type { return log(el); } ;
  
  data_array test_array_hand_logged = test_array;
  for(int i=0; i< test_array.get_dims(0); i++){
    for(int j=0; j< test_array.get_dims(1); j++){
      test_array_hand_logged.set_element(i, j, log(test_array_hand_logged.get_element(i, j)));
    }
  }
  test_array.apply(log_function);
  if(compare_2d(test_array, test_array_hand_logged)){
  
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Apply function error", 1);
  }
  if(err == TEST_PASSED) test_bed->report_info("Basic tests OK");
  return err;
}

int test_entity_data_array::three_d_and_shift(){
/** Test 3-d resizer and shift*/

  int err= TEST_PASSED;
  int tmp_err;
  //And now a 3-d version
  test_array = data_array(10, 10, 10);
  for(size_t i=0; i<test_array.get_dims(0); i++){
    for(size_t j =0; j<test_array.get_dims(1); j++){
      for(size_t k =0; k<test_array.get_dims(2); k++){
        tmp_err=test_array.set_element(i, j, k, (i+1)*(2*j+1)*(4*k+1));
        if(tmp_err) err |= TEST_ASSERT_FAIL;
      }
    }
  }
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

  //Now test the shift function

  test_array.shift(1, 3);
  test_array.shift(1, -3);

  if(compare_3d(old_array, test_array, true)){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }

  test_array.shift(2, 3);
  test_array.shift(2, -3);

  if(compare_3d(old_array, test_array, true)){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }

  test_array.shift(0, 2);

  if(!compare_3d(old_array, test_array, true)){
  //If they still compare equal we has problem
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, no shift applied", 1);
  }

  test_array.shift(0, -2);
  
  if(compare_3d(old_array, test_array, true)){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Shift error, wrong values read", 1);
  }
  return err;

}

int test_entity_data_array::technical_tests(){
/** Check things like copy constructors
*
* Expects a 2-d array and wont be any use if els and axes not set
*/

  test_bed->report_info("Checking technical aspects ", 2);

  test_array = data_array(10, 10);

  int err = TEST_PASSED;
  data_array dat = test_array;
  err |= compare_2d(test_array, dat);
  
  try{
    std::vector<data_array> my_vec;
    my_vec.push_back(test_array);
    my_vec.push_back(test_array);
    my_vec.resize(100);
    my_vec.resize(10);
    err|= compare_2d(my_vec[0], test_array);
  }catch(const std::exception& e){
    //Swallow and continue if possible,
    std::string message = e.what();
    test_bed->report_info("Exception message " +message, 1);
    err |= TEST_ASSERT_FAIL;
  }
  return err;
}

int test_entity_data_array::io_tests(){

  test_bed->report_info("Checking file io", 1);
  int err = TEST_PASSED;
  bool err2=false;
  std::string filename = "./files/test_file.dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    err2=test_array.write_to_file(file);
  }
  file.close();
  if(err2) test_bed->report_info("Error writing testfile", 1);
  
  data_array new_array = test_array;
  file.open(filename.c_str(),std::ios::in|std::ios::binary);
  err2=new_array.read_from_file(file);
  if(err2) test_bed->report_info("Error reading testfile", 1);
  file.close();
  
  if(err2) err|=TEST_ASSERT_FAIL;
  if(compare_3d(test_array, new_array)) err |= TEST_WRONG_RESULT;

  return err;
}

//----------------------------------------------------------------

test_entity_get_and_fft::test_entity_get_and_fft(){
  name = "read and FFT";

}

test_entity_get_and_fft::~test_entity_get_and_fft(){
 
  if(test_rdr) delete test_rdr;
  
}

int test_entity_get_and_fft::run(){
/** \brief Checks full read and fft procedure
*
*Reads a test sdf file, stores into data array and runs fft. Test data should be a sine curve with one major frequency which is then checked. Note frequency is hard coded to match that produced by ./files/sin.deck
*/

  int err = TEST_PASSED;
  
  char block_id[ID_SIZE]= "ex";
  test_rdr = new reader("./files/sin", block_id);

  test_dat = data_array(10, 10);

  err |=one_d();

  //strcpy(block_id, "ay");
  strcpy(block_id, "ax");

  if(test_rdr) delete test_rdr;
  test_rdr = new reader("./files/sinAcc", block_id);
  err|= two_d();

  test_bed->report_err(err);
  return err;

}

int test_entity_get_and_fft::one_d(){
  int err = TEST_PASSED;

  int tim_in[3], space_in[2];
  tim_in[0]=0;
  tim_in[1]=1;
  tim_in[2]=0;
  space_in[0]=0;

  int n_tims = std::max(tim_in[1]-tim_in[0], 1);

  size_t n_dims;
  std::vector<size_t> dims;
  test_rdr->read_dims(n_dims, dims);

  space_in[1]=dims[0];

  if(n_dims !=1){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Array dims wrong", 1);
    return err;
  }

  test_dat = data_array(dims[0], n_tims);
  test_dat_fft = data_array(dims[0], n_tims);
  if(!test_dat.is_good()||!test_dat_fft.is_good()){
    err|=TEST_ASSERT_FAIL;
    return err;
  }
  
  test_rdr->read_data(test_dat, tim_in, space_in);

  {
    bool tmp_err = test_dat.fft_me(test_dat_fft);
    if(tmp_err) err|=TEST_ASSERT_FAIL;
    if(test_dat_fft.check_ids(test_dat)) err |= TEST_WRONG_RESULT;
    if(err == TEST_PASSED) test_bed->report_info("1D read and FFT reports no error", 1);

    //Get primary frequency
    int max_index = 0;
    my_type max_val = 0;
    std::vector<size_t> max_pos;
    my_type expected_max = 1.2566371e-4;
    bool both_freqs_correct = true;

    max_val = test_dat_fft.maxval(max_pos);
    if(max_pos.size() <1) err |=TEST_WRONG_RESULT;
    max_index = max_pos[0];
    if(std::abs(std::abs(test_dat_fft.get_axis_element(0,max_index)) - expected_max) > PRECISION){
      err|= TEST_WRONG_RESULT;
      test_bed->report_info("Max freq is "+mk_str(test_dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+")", 1);
      both_freqs_correct = false;
    }

    max_val = test_dat_fft.maxval(max_pos, max_index+1);
    if(max_pos.size() <1) err |=TEST_WRONG_RESULT;
    max_index = max_pos[0];
    if(std::abs(std::abs(test_dat_fft.get_axis_element(0,max_index)) - expected_max) > PRECISION){
      err|= TEST_WRONG_RESULT;
      test_bed->report_info("Max freq is "+mk_str(test_dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+")", 1);
      both_freqs_correct = false;
    }
    
    if(both_freqs_correct) test_bed->report_info("FFT Frequency correct!", 1);
  }
  //Now size down by 1 element and redo. This checks odd and even total sizes
  {

    test_dat.resize(0, dims[0]-1, true);
    test_dat_fft.resize(0, dims[0]-1, true);
    bool tmp_err = test_dat.fft_me(test_dat_fft);

    if(tmp_err) err|=TEST_ASSERT_FAIL;
    if(test_dat_fft.check_ids(test_dat)) err |= TEST_WRONG_RESULT;
    if(err == TEST_PASSED) test_bed->report_info("1D read and FFT reports no error", 1);

    //Get primary frequency
    int max_index = 0;
    my_type max_val = 0;
    std::vector<size_t> max_pos;
    my_type expected_max = 1.2566371e-4;
    bool both_freqs_correct = true;

    max_val = test_dat_fft.maxval(max_pos);
    if(max_pos.size() <1) err |=TEST_WRONG_RESULT;
    max_index = max_pos[0];
    int sgn = test_dat_fft.get_axis_element(0,max_index)/std::abs(test_dat_fft.get_axis_element(0,max_index));
    if( !((test_dat_fft.get_axis_element(0,max_index+1) > sgn*expected_max) && (test_dat_fft.get_axis_element(0,max_index-1) < sgn*expected_max))){
      err|= TEST_WRONG_RESULT;
      test_bed->report_info("Max freq is "+mk_str(test_dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+")", 1);
      both_freqs_correct = false;
    }
  
    max_val = test_dat_fft.maxval(max_pos, max_index+1);
    if(max_pos.size() <1) err |=TEST_WRONG_RESULT;
    max_index = max_pos[0];
    sgn = test_dat_fft.get_axis_element(0,max_index)/std::abs(test_dat_fft.get_axis_element(0,max_index));
    if(!( (test_dat_fft.get_axis_element(0,max_index+1) > sgn*expected_max) && (test_dat_fft.get_axis_element(0,max_index-1) < sgn*expected_max))){
      err|= TEST_WRONG_RESULT;
      test_bed->report_info("Max freq is "+mk_str(test_dat_fft.get_axis_element(0,max_index))+" ("+mk_str(max_index)+")", 1);
      both_freqs_correct = false;
    }
    
    if(both_freqs_correct) test_bed->report_info("FFT Frequency correct!", 1);
  }
  
  return err;

}

int test_entity_get_and_fft::two_d(){
  int err = TEST_PASSED;

  int tim_in[3], space_in[2];
  tim_in[0]=0;
  tim_in[1]=3;
  tim_in[2]=100;
  space_in[0]=0;

  int n_tims = tim_in[2];//std::max(tim_in[1]-tim_in[0], 1);

  int space_size;
  size_t n_dims;
  std::vector<size_t> dims;
  test_rdr->read_dims(n_dims, dims);

  space_in[1]=dims[0];
  space_size = space_in[1]-space_in[0];
  if(n_dims !=1){
    err |= TEST_WRONG_RESULT;
    test_bed->report_info("Array dims wrong", 1);
    test_bed->report_err(err);

    return err;
    //nothing more worth doing right now...
  }

  test_dat = data_array(space_size, n_tims);
  test_dat_fft = data_array(space_size, n_tims);
  if(!test_dat.is_good()||!test_dat_fft.is_good()){
    err|=TEST_ASSERT_FAIL;
    return err;
  }
  
  test_rdr->read_data(test_dat, tim_in, space_in);

  bool tmp_err = test_dat.fft_me(test_dat_fft);
  if(tmp_err) err|=TEST_ASSERT_FAIL;
  if(test_dat_fft.check_ids(test_dat)) err |= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("2D read and FFT reports no error", 1);
  
  int max_index = 0;
  my_type max_val = 0;
  std::vector<size_t> max_pos;
  my_type expected_max = 1.2566371e-4;
  bool both_freqs_correct = true;

  max_val = test_dat_fft.maxval(max_pos);

  if(max_pos.size() <2) err |=TEST_WRONG_RESULT;
  max_index = max_pos[0];
//  for(int i=0; i<max_pos.size(); i++) std::cout<<max_pos[i]<<" ";
  
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
  
  filename = test_rdr->file_prefix+"FFT2DTest_"+block +"_"+time_str+"_"+mk_str(test_dat_fft.space[0])+"_"+mk_str(test_dat_fft.space[1]) + ".dat";
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
/** \todo Add setup and teardown so these are alloced before run but not on construction*/
  name = "basic maths helpers";
  size = 256;
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
test_entity_basic_maths::~test_entity_basic_maths(){

  free(data_square);
  free(data_positive);
  free(data_tmp);
  free(axis);
  free(d_axis);
  free(axisf);
}

int test_entity_basic_maths::run(){
  int err=TEST_PASSED;

  my_type target;
  int whe;

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
    //test_bed->report_info(mk_str(target)+" "+mk_str(axisf[whe]), 2);
  }
  if(err == TEST_PASSED) test_bed->report_info("Where OK", 1);

  calc_type res = integrator(data_square, size, d_axis);
  if(res!= 0.0) err |= TEST_WRONG_RESULT;
  res = integrator(data_positive, size, d_axis);
  if(std::abs(res - (pow((calc_type)(size-1), 2)/20.0)) > res*PRECISION) err |= TEST_WRONG_RESULT;
  //test it's correct to within some finite precision, defined in header
  if(err == TEST_PASSED) test_bed->report_info("Integrator OK", 1);

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
    //std::cout<<data_positive[i]<<" "<<data_tmp[i]<<'\n';
    total += std::abs(data_positive[i] - data_tmp[i]);
  }
  //Smooth on straight line with odd width should do nothing except at ends...
  if(std::abs(total) > PRECISION) err |=TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Boxcar smooth OK", 1);


  {
    calc_type axis[4]={0.0,1.0,2.0,3.0}, vals[4]={0.0,1.0,0.0,1.0}, target, interp;

    target = 0.25;
    interp = interpolate(axis, vals, target, 1);
    if(std::abs(interp - 0.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 0.75;
    interp = interpolate(axis, vals, target, 1);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    //Nearest value interpol
    target = 0.5;
    interp = interpolate(axis, vals, target, 2);
    if(std::abs(interp - 0.5)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.0;
    interp = interpolate(axis, vals, target, 2);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.5;
    interp = interpolate(axis+1, vals+1, target, 2);
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

  test_bed->report_err(err);
  return err;

}

test_entity_extern_maths::test_entity_extern_maths(){

  name = "external maths routines";

}
test_entity_extern_maths::~test_entity_extern_maths(){

}

int test_entity_extern_maths::run(){

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
  //test_bed->report_info(mk_str(std::abs(bess + bess1 - 2.0*(double)index*bess2/arg)), 1);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;

  index=7;
  arg=-20.98;
  bess = boost::math::cyl_bessel_j(index-1, arg);
  bess1 = boost::math::cyl_bessel_j(index+1, arg);
  bess2 = boost::math::cyl_bessel_j(index, arg);
  //test_bed->report_info(mk_str(std::abs(bess + bess1 - 2.0*(double)index*bess2/arg)), 1);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;


  index =0;
  arg =1.0;
  bess = boost::math::cyl_bessel_j(index, arg);

  if(std::abs(bess - 0.7651976865579665514497) >PRECISION) err|= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Bessel functions OK", 1);

  test_bed->report_err(err);

  return err;

}
//----------------------------------------------------------------

test_entity_plasma::test_entity_plasma(){
  name = "plasma";
  plas = new plasma("./files/test");

}
test_entity_plasma::~test_entity_plasma(){

  delete plas;
}

int test_entity_plasma::run(){
/** \brief Test resonant frequencies and refractive indices
*
*Checks the resonant frequencies obey the equations used to derive them. Checks the dispersion roots for Whistlers match those found using high-density approx. Checks plasma O and X mode dispersion too. Note, first call for issues with these tests is to check returned mu.err on failing tests
*/

  int err=TEST_PASSED;
  
  err |= analytic_dispersion();
  err |= resonant_freq();
  err |= high_density();
  err |= other_modes();
  err |= phi_dom();

  test_bed->report_err(err);
  return err;

}

int test_entity_plasma::analytic_dispersion(){
/** \brief Check analytic dispersion relations
*
*Checks the analytic relations, both ways and including derivatives \todo Check theta
*/

  int err=TEST_PASSED;

  calc_type k, om, om_new, d_om, d_k;
  size_t n_tests = 5;
  //First whistlers
  //At very large k we should get om_ce and at 0, we get 0
  k=100;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  if(std::abs(om /plas->get_omega_ref("ce") -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  k = 0;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  if(std::abs(om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  //Symmetry
  calc_type om_2;
  k=0.0001;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  k=-k;
  om_2 = plas->get_dispersion(k, WAVE_WHISTLER);
  
  if(std::abs(om - om_2) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  
  //If we pick a few e.g.s and do them both ways we should get same result back. Start with omega for simplicity
  for(size_t i=1; i<= n_tests; i++){
    om = (float)i *plas->get_omega_ref("ce")/(float) (n_tests +1);
    k = plas->get_dispersion(om, WAVE_WHISTLER, true);
    om_new = plas->get_dispersion(k, WAVE_WHISTLER);
    if(std::abs(om /om_new -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Derivs should be 0 at 0 and large omega.

  k=100;
  d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  k = 0;
  d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  
  //Finally derivs at our samples should be inverses
  for(size_t i=1; i<= n_tests; i++){
    om = (float)i *plas->get_omega_ref("ce")/(float) (n_tests +1);
    k = plas->get_dispersion(om, WAVE_WHISTLER, true);
    
    d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
    d_k = plas->get_dispersion(om, WAVE_WHISTLER, true, true);
    if(std::abs(d_om*d_k -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Finally we check how it handles out of range
  try{
    k = plas->get_dispersion(plas->get_omega_ref("ce")*1.5, WAVE_WHISTLER, true);
  }catch(...){
    err |= TEST_ASSERT_FAIL;
  }

//  std::cout<<"Err: "<<err<<'\n';
  //Now EM
  //At 0, we get om_pe
  k = 0;
  om = plas->get_dispersion(k, WAVE_O);
  if(std::abs(om - plas->get_omega_ref("pe")) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  
  //If we pick a few e.g.s and do them both ways we should get same result back. Go from say om_pe to 3om_pe
  for(size_t i=1; i<= n_tests; i++){
    om = plas->get_omega_ref("pe")*(1.0+ + 2.0*(float)i/(float) (n_tests +1));
    k = plas->get_dispersion(om, WAVE_O, true);
    om_new = plas->get_dispersion(k, WAVE_O);
    if(std::abs(om /om_new -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Derivs should be 0 at 0 and ~1/c at very large omega
  k = 0;
  d_om = plas->get_dispersion(k, WAVE_O, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  om = 100* plas->get_omega_ref("pe");
  d_k = plas->get_dispersion(om, WAVE_O, true, true);
  if(std::abs(d_k*v0 - 1.0) > LOW_PRECISION) err |=TEST_WRONG_RESULT;
  std::cout<<"Err: "<<err<<'\n';

  //Finally derivs at our samples should be inverses
  for(size_t i=1; i<= n_tests; i++){
    om = plas->get_omega_ref("pe")*(1.0+ + 2.0*(float)i/(float) (n_tests +1));
    k = plas->get_dispersion(om, WAVE_O, true);
    
    d_om = plas->get_dispersion(k, WAVE_O, false, true);
    d_k = plas->get_dispersion(om, WAVE_O, true, true);
    if(std::abs(d_om*d_k -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  //Finally we check how it handles out of range
  try{
    k = plas->get_dispersion(plas->get_omega_ref("pe")*0.8, WAVE_O, true);
  }catch(...){
    err |= TEST_ASSERT_FAIL;
  }

  if(err == TEST_PASSED) test_bed->report_info("Analytic dispersion OK", 1);
  else  test_bed->report_info("Error in analytic dispersion", 1);

  return err;
}

int test_entity_plasma::resonant_freq(){
/** \brief Check resonant frequency solver
*
*Checks the returned resonant frequency obeys equations used to derive it by solving both for mu.
*/

  int err=TEST_PASSED;
  int n_tests = 5;
  std::vector<calc_type> results;
  mu_dmudom my_mu;
  int err_count = 0;
  calc_type x, v_par, n, om_ce_local, om_pe_local;
  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");


  calc_type cos_theta, mu_tmp1, mu_tmp2, tmp_omega_n=0;
  calc_type gamma, gamma2;

  test_bed->report_info("Testing resonant frequency solver", 1);


  for(int ii=0; ii<n_tests; ii++){
    v_par = (0.01 + 0.5*(float)ii/ (float)(n_tests+1))* v0;

    for(int j=0; j< n_tests; j++){
      x = 4.0* (float) j / (float)(n_tests+1);
      cos_theta = std::cos(std::atan(x));

      for(int k=0; k< n_tests; k++){
        n = -n_tests/2 + k*n_tests/2;
        
        gamma2 = 1.0/( 1.0 - std::pow(v_par/v0, 2));
        
        gamma = std::sqrt(gamma2);
        
        results = plas->get_resonant_omega(x, v_par, n);
        /**Now check each element of the resonant frequency solution set satisfies Stix 2.45 and the resonance condition together*/
        for(int i=0; i<(int)results.size(); ++i){
          //test_bed->report_info("Freq is "+mk_str(results[i], true)+" = "+mk_str(results[i]/my_const.omega_ce, true)+" om_ce", 2);
          
          mu_tmp1 = std::pow(v0 * (gamma*results[i] - n*om_ce_local)/(gamma*results[i] * v_par *cos_theta), 2);
          mu_tmp2 = (1.0 - (std::pow(om_pe_local,2)/(results[i]*(results[i] + om_ce_local*cos_theta))));

          if(std::abs((mu_tmp1 - mu_tmp2)/mu_tmp1) > NUM_PRECISION){
            err|=TEST_WRONG_RESULT;
            test_bed->report_info("refractive index mismatch of "+mk_str((mu_tmp1-mu_tmp2)/mu_tmp1), 2);
          }
        
          //Also check there is a valid full mu solution
          my_mu = plas->get_high_dens_phi_mu_om(results[i], std::atan(x), 0.0, 0.0, tmp_omega_n);
          if(my_mu.err){
            err|=TEST_WRONG_RESULT;
            test_bed->report_info("No full mu solution for resonant frequency", 2);
            err_count++;
          }
        }
      }
    }
  }
  test_bed->report_info("Mu error count: "+mk_str(err_count), 2);
  return err;

}

int test_entity_plasma::high_density(){
/** \brief Tests high density approximation for dispersion relations
*
*Test if the mu found by get_root and get_phi_mu_om matches the high density whistler in high dens regime
*/

  int err=TEST_PASSED;

  calc_type om_ce_local, om_pe_local;
  calc_type mu_tmp2;

  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");

  size_t n_tests = 10;
  calc_type tmp_omega=0.0, tmp_theta=pi/(calc_type)(n_tests), tmp_omega_n=0.0;
  mu_dmudom my_mu;
  mu my_mu_all;
  mu_dmudom my_mu_dens;
  int err_cnt=0;
  test_bed->report_info("Testing whistler high density approx.", 1);

  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_ce_local)/(calc_type)(n_tests + 1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_dens = plas->get_high_dens_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /** my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
      test_bed->report_info("Mismatch in high density approx or dispersion solver at "+mk_str(tmp_omega/std::abs(om_ce_local))+" "+mk_str(tmp_theta), 1);
      test_bed->report_info("Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu, true), 2);
    }
    /** my_mu_dens should EXACTLY equal Stix 2.45 without the 1.0 term*/
    mu_tmp2 = sqrt( - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));

    if(std::abs(my_mu_dens.mu-mu_tmp2)/my_mu_dens.mu > NUM_PRECISION){
      err_cnt++;
      test_bed->report_info("    Mismatch in alternate dispersion solver at "+mk_str(tmp_omega/std::abs(om_ce_local))+" "+mk_str(tmp_theta), 1);
      test_bed->report_info("    Mu "+mk_str(my_mu_dens.mu)+" difference "+mk_str(my_mu_dens.mu - mu_tmp2)+" relative error "+mk_str((my_mu_dens.mu-mu_tmp2)/my_mu_dens.mu), 2);
    }

    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }
  
  tmp_omega = 0.6*std::abs(om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_theta += pi/(calc_type)(n_tests);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /* my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
    
      test_bed->report_info("Mismatch in high density approx or dispersion solver at "+mk_str(tmp_omega/std::abs(om_ce_local))+" "+mk_str(tmp_theta), 1);
      test_bed->report_info("Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu), 2);
    }
     //my_mu_all.mu and my_mu.mu should be exactly equal:
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }

 
  }
  if(err_cnt> 0){
    test_bed->report_info("Total "+mk_str(err_cnt)+" out of "+mk_str(2*(int)n_tests)+" issues in high density approx or dispersion solver at precision: "+mk_str(LOW_PRECISION), 1);
    //err|=TEST_WRONG_RESULT;
    //Make these a warning not an error because we expect them sometimes
  }

  return err;
}

int test_entity_plasma::other_modes(){
/** \brief Test dispersion solver with other wave modes */

  int err=TEST_PASSED;

  calc_type om_ce_local, om_pe_local;
  calc_type mu_tmp2;

  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");


  test_bed->report_info("Testing dispersion solver for plasma O mode", 1);
  size_t n_tests = 10;
  mu_dmudom my_mu;
  mu my_mu_all;

  /**Try plasma wave modes in solvers, perpendicular propagation*/
  calc_type tmp_omega = om_pe_local, tmp_omega_n=0;
  calc_type tmp_theta = pi/2.0;
  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_pe_local)/(calc_type)(n_tests + 1);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);

    mu_tmp2 = std::sqrt(std::pow(tmp_omega, 2) - std::pow(om_pe_local, 2))/tmp_omega;
    
    if(std::abs(my_mu_all.mu-mu_tmp2)/my_mu_all.mu > LOW_PRECISION){
      
      test_bed->report_info("Error in approx or dispersion solver for plasma wave at "+mk_str(tmp_omega/std::abs(om_pe_local))+" "+mk_str(tmp_theta), 1);
      
      test_bed->report_info("Mu "+mk_str(my_mu_all.mu)+" difference "+mk_str(my_mu_all.mu - mu_tmp2)+" relative error "+mk_str((my_mu_all.mu-mu_tmp2)/my_mu_all.mu), 2);
    }
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;

    }
    
  }
  /**Try left hand X mode too*/
  test_bed->report_info("Testing dispersion solver for plasma X mode", 1);

  calc_type omega_UH = std::sqrt(om_pe_local*om_pe_local + om_ce_local*om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_pe_local)/(calc_type)(n_tests + 1);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta, false);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n, false);
    
    mu_tmp2 = std::sqrt(1.0 - std::pow(om_pe_local/tmp_omega, 2)*(std::pow(tmp_omega, 2) - std::pow(om_pe_local, 2))/(std::pow(tmp_omega, 2) - std::pow(omega_UH, 2)));
    if(std::abs(my_mu_all.mu-mu_tmp2)/my_mu_all.mu > LOW_PRECISION){
      
      test_bed->report_info("Error in approx or dispersion solver for plasma wave at "+mk_str(tmp_omega/std::abs(om_pe_local))+" "+mk_str(tmp_theta), 1);
      
      test_bed->report_info("Mu "+mk_str(my_mu_all.mu)+" difference "+mk_str(my_mu_all.mu - mu_tmp2)+" relative error "+mk_str((my_mu_all.mu-mu_tmp2)/my_mu_all.mu), 2);
    }
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;

    }
    
  }

  return err;

}

int test_entity_plasma::phi_dom(){
/** \brief Test other plasma returns
*
* Checks the values of mu.dom, mu.dmudtheta and phi against special cases. \todo Write this Phi is pretty utestable really...
*/

  int err=TEST_PASSED;
  size_t n_tests = 10;

  calc_type mu_tmp1, mu_tmp2, om_ce_local, om_pe_local;
  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");

  calc_type d_omega = std::abs(om_ce_local)/1e8;
  calc_type d_theta = pi/(calc_type)(n_tests)/1e7;
  //Derivative step size.

  mu_dmudom my_mu, my_mu_p;
  mu my_mu_all, my_mu_all_p;

  calc_type tmp_omega = 0.0, tmp_theta=pi/(calc_type)(n_tests), tmp_omega_n=0.0;

  test_bed->report_info("Testing dmu/domega", 1);

  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_ce_local)/(calc_type)(n_tests + 1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu_p = plas->get_phi_mu_om(tmp_omega+d_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all_p = plas->get_root(0.0, tmp_omega+d_omega, tmp_theta);

    /** Approx numerical derivative*/
    mu_tmp1 = (my_mu.mu - my_mu_p.mu)/d_omega;
    mu_tmp2 = (my_mu_all.mu - my_mu_all_p.mu)/d_omega;
    if(std::abs(std::abs(mu_tmp1/my_mu.dmudom) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_phi_mu_om", 2);
    }
    if(std::abs(std::abs(mu_tmp2/my_mu_all.dmudom) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_root", 2);
    }

    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.dmudom-my_mu.dmudom) > PRECISION){
      test_bed->report_info("Inconsistent derivative between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }
  
  test_bed->report_info("Testing dmu/dtheta", 1);
  
  tmp_omega = 0.6*std::abs(om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_theta += pi/(calc_type)(n_tests);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu_p = plas->get_phi_mu_om(tmp_omega, tmp_theta+d_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all_p = plas->get_root(0.0, tmp_omega, tmp_theta+d_theta);

    /** Approx numerical derivative. Manually fix signs*/
    mu_tmp1 = -(my_mu.mu - my_mu_p.mu)/d_theta;

    mu_tmp2 = -(my_mu_all.mu - my_mu_all_p.mu)/d_theta;
    
    if(std::abs(std::abs(mu_tmp1 /my_mu.dmudtheta) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_phi_mu_om at omega = "+mk_str(tmp_omega, true) +" and phi = "+mk_str(tmp_theta, true), 2);

    }
    if(std::abs(std::abs(mu_tmp2/my_mu_all.dmudtheta) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_root at omega = "+mk_str(tmp_omega, true) +" and phi = "+mk_str(tmp_theta, true), 2);

    }
    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.dmudtheta-my_mu.dmudtheta) > PRECISION){
      test_bed->report_info("Inconsistent derivative between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }

  return err;

}

/*
For getting the resonant frequency to consider errors will matter, but the noise in spectra can be expected to be a similar sort of error
*
The reason for using the better dispersion solver is a) to avoid any numerical derivatives in general and b) because some of the calculation depends on derivs of the refractive index

*So basically I am using the easy approx where a) it’s a genuine nightmare to do better (10-14th order polynomial) and b) I expect other errors to be similar or more important and c) I really hope it’s linear or polynomial in error*/
//----------------------------------------------------------------

test_entity_spectrum::test_entity_spectrum(){

  name = "spectrum checks";
  file_prefix = "./files/";
  
}
test_entity_spectrum::~test_entity_spectrum(){

  if(test_contr) delete test_contr;

}

int test_entity_spectrum::run(){
/** \brief Test spectrum extraction
*
* This should test the dispersion relation approximations are OK (plain and vg). Check test spectrum makes sense. Test extraction of a spectrum from data. Note data does not come from files, but from a test file which is already written as a data array using ./files/generate_fftd.pro which makes FFT_data.dat with the FFTd data and spectrum.dat with a derived spectrum to check against.
*/

  int err = TEST_PASSED;

  //Use a different deck.status file...
  deck_constants const_tmp = my_const;
  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
  share_consts();

  err|= setup();
  if(!test_bed->check_for_abort(err)) err|= basic_tests1();
  if(!test_bed->check_for_abort(err)) err|= basic_tests2();
  if(!test_bed->check_for_abort(err)) err|= albertGs_tests();
  
  test_bed->report_err(err);

  my_const = const_tmp;
  share_consts();

  return err;

}

int test_entity_spectrum::setup(){
/** \brief Setup to test spectrum
*
* Note strictly this is the test of data array constructor taking a filename too.
*/

  int err = TEST_PASSED;

  test_dat_fft = data_array(file_prefix + "FFT_data.dat", true);
  if(!test_dat_fft.is_good()){
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }
  test_contr = new controller(file_prefix);

  return err;
}

int test_entity_spectrum::basic_tests1(){
/** \brief Basic tests of spectrum
*
*Test test_spectrum and angle generation \todo The angles are integrating to 0.5 not 1. Which do we want????
*/
  int err = TEST_PASSED;

  std::fstream outfile;
  size_t len=0;
  my_type total_error =0.0;
  my_type * d_angle, * angle_data;

  test_contr->add_spectrum(test_dat_fft.get_dims(0), DEFAULT_N_ANG, true);

  /** Check this test spectrum makes sense \todo HOW????*/

  test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_DELTA);
  //Check angle distrib integrates to 1 for each case
  //NOTE we can only do this if MIN_ANG is either 0 or is - MAX_ANG. otherwise we're into erf and bunk
  bool is_symmetric=false, is_zero = false;
  if(std::abs(ANG_MIN + ANG_MAX) < PRECISION) is_symmetric = true;
  if(std::abs(ANG_MIN) < PRECISION) is_zero = true;
  
  if(is_symmetric || is_zero){
    len = test_contr->get_current_spectrum()->get_g_dims(1);
    d_angle = (my_type *) calloc(DEFAULT_N_ANG, sizeof(my_type));
    for(size_t i=0; i<DEFAULT_N_ANG-1; ++i){
      d_angle[i] = std::abs(test_contr->get_current_spectrum()->get_ang_axis_element(i) - test_contr->get_current_spectrum()->get_ang_axis_element(i+1));
    }
    angle_data = (my_type *) malloc(len*sizeof(my_type));
    for(size_t i=0; i<len; i++){
      *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);
    }
    
    total_error = integrator(angle_data, len, d_angle);
    test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_GAUSS);
    for(size_t i=0; i<len; i++) *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);

    total_error += integrator(angle_data, len, d_angle);

    test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_ISO);
    for(size_t i=0; i<len; i++) *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);
    total_error += integrator(angle_data, len, d_angle);
    
    my_type expected = is_zero ? 2.0 : 3.0;
    //Iso always integrates to 1. Gaussian and delta are always symmetric
    if(std::abs(total_error - expected)/3.0 > NUM_PRECISION){
    
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("Error in angular distribution integrals, value " + mk_str(total_error, true));
    }
    if(angle_data) free(angle_data);
    if(d_angle) free(d_angle);
  }else{
    test_bed->report_info("Cannot test assymmetric spectrum");
  }
  outfile.open("spect_testy.dat", std::ios::binary|std::ios::out);
  test_contr->get_current_spectrum()->write_to_file(outfile);
  outfile.close();
  if(err == TEST_PASSED) test_bed->report_info("Test spectrum OK");
  
  /** Now make the real spectrum from data and check the result matches the plain text test file*/
  return err;
}
int test_entity_spectrum::basic_tests2(){
/**\brief Test spectrum extraction
*
*Compare extracted spectrum from FFT'd data file to test file
*/
  int err = TEST_PASSED;

  std::fstream outfile, infile;
  size_t len=0;
  my_type total_error =0.0;

  test_contr->get_current_spectrum()->generate_spectrum(test_dat_fft ,10, FUNCTION_GAUSS);

  test_spect = data_array(file_prefix + "spectrum.dat", true);
  if(test_spect.is_good()){
    //We ignore frequencies below say 0.05 om_ce
    my_type * ax = test_spect.get_axis(0, len);
    int min_ind = 0;
    if(ax) min_ind = where(ax+len/2, len/2, 17588.200*0.05);
    /**Hard code min freq to match the IDL file with test data generation...*/
    
    total_error = 0.0;
    for(size_t i=0; i< len/2 - min_ind; i++){
      total_error += std::abs(test_contr->get_current_spectrum()->get_B_element(i)-test_spect.get_element(i));
    }
    for(size_t i=len/2 + min_ind; i< len; i++){
      total_error += std::abs(test_contr->get_current_spectrum()->get_B_element(i)-test_spect.get_element(i));

    }
    if(total_error > LOW_PRECISION){
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("Mismatch between generated spectrum and test spectrum of "+mk_str(total_error));
    }
    /* Preserve the spectrum*/
    outfile.open("spect_out.dat", std::ios::binary|std::ios::out);
    test_contr->get_current_spectrum()->write_to_file(outfile);
    outfile.close();
  }else{
    err |= TEST_ASSERT_FAIL;
  }
  if(err == TEST_PASSED) test_bed->report_info("Generate spectrum OK");
  
  data_array old_B = test_contr->get_current_spectrum()->copy_out_B();
  //Now dump to file and read back in and compare
  bool err2 = test_contr->add_spectrum("spect_out.dat");
  if(err2) err |= TEST_ASSERT_FAIL;
  else{
    data_array new_B = test_contr->get_current_spectrum()->copy_out_B();
    if(!old_B.is_good() || !new_B.is_good() || compare_2d(old_B, new_B)){
      test_bed->report_info("Error or Mismatch in read", 0);
      err|= TEST_WRONG_RESULT;
    }
  }
  return err;

}

int test_entity_spectrum::albertGs_tests(){
/** \brief Tests of the Albert G functions in spectrum. Also tests the normalisations on the way. NOTE: since we're comparing the values of an analytic function with a numerical integral, we can get mismatches at the cutoffs. More points should help this. If that doesn't there may be something wrong.
*
*
*/
  int err = TEST_PASSED;

  calc_type om_ce_local, om_pe_local, G1, G2, G1_analytic, G2_analytic;
  om_ce_local = test_contr->get_plasma().get_omega_ref("ce");
  om_pe_local = test_contr->get_plasma().get_omega_ref("pe");

  calc_type mass_ratio = 1836.2;

  size_t n_tests = 10;
  calc_type tmp_omega=0.0, tmp_x;

  test_contr->add_spectrum(5000, DEFAULT_N_ANG, true);

  if(!test_contr->get_current_spectrum()->is_good()){
    my_print("Spectrum in invalid state. Aborting", mpi_info.rank);
    err |=TEST_ASSERT_FAIL;
    err |=TEST_FATAL_ERR;
    return err;
  }

  test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_GAUSS);
  
  my_type om_min, om_max, x_min, x_max, om_peak;
  om_min = 12000.0;
  om_max = 16500.0;
  //make sure this is lower than the test spectrum axis range
  x_min = 0.0;
  x_max = 3.0;
  
  test_contr->get_current_spectrum()->truncate_om(om_min, om_max);
  test_contr->get_current_spectrum()->truncate_x(x_min, x_max);
  om_peak = test_contr->get_current_spectrum()->get_peak_omega();
  //Now we have a test spectrum. Need to know what its normalisations should be. And what the Albert functions should resolve to.
  
  my_type width=0.1*om_peak;
  
  for(size_t i=0; i< n_tests;i++){
    tmp_omega = std::abs(om_ce_local)/10.0 + 89.0/100.0 * om_max * (1.0 - exp(-i));
    //Cover range from small to just below om_ce...
    G1 = test_contr->get_current_spectrum()->get_G1(tmp_omega);

    //Analytic calculations for truncated Gaussians, see Albert 2005
    
    if(tmp_omega > om_min && tmp_omega < om_max){
      G1_analytic = 2.0 / std::sqrt(pi) * std::exp( - std::pow((tmp_omega - om_peak)/width, 2));
      G1_analytic /= (boost::math::erf((om_max - om_peak)/width) +boost::math::erf((om_peak - om_min)/width));
      G1_analytic /=width;
    }else{
      G1_analytic = 0.0;
    }

    if( (G1 != 0.0 && std::abs(G1-G1_analytic)/(G1) > LOW_PRECISION)|| (G1 == 0.0 && G1_analytic != 0.0)){
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("G1 does not match analytic calc, relative error = "+mk_str((std::abs(G1/G1_analytic)-1.0)*100, true)+" at "+mk_str(tmp_omega, true), mpi_info.rank);
    }
  }


  tmp_omega = 0.6 * std::abs(om_ce_local);
  for(size_t i=0; i< n_tests;i++){

    tmp_x = ANG_MIN + i * (ANG_MAX - ANG_MIN)/(n_tests-1);
    
    G2 = test_contr->get_current_spectrum()->get_G2(tmp_omega, tmp_x);
    
  //  std::cout<<"  G2 "<<G2<<std::endl;
    G2_analytic = std::pow((( mass_ratio / (1.0 + mass_ratio))*om_ce_local*om_ce_local/om_pe_local/om_pe_local), 1.5);
    G2_analytic *= exp(- (tmp_x*tmp_x)/std::pow(SPECTRUM_ANG_STDDEV, 2));
    
    }


  std::fstream outfile;
  outfile.open("spect_truncated.dat", std::ios::out|std::ios::binary);
  test_contr->get_current_spectrum()->write_to_file(outfile);
  outfile.close();


  return err;
}
//----------------------------------------------------------------

test_entity_levelone::test_entity_levelone(){

  name = "level-one derivation";
  
  test_contr=nullptr;
  my_reader=nullptr;
  
}
test_entity_levelone::~test_entity_levelone(){

  delete test_contr;
  delete my_reader;
}

int test_entity_levelone::run(){
/** \brief Test entire level-1 data extraction
*
**/

  int err = TEST_PASSED;

  //Use a different deck.status file...
//  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
//  share_consts();

  strcpy(block_id, "ay");

  file_prefix = "./files/l1/l1";
  space_in[0] = 0;
  space_in[1] = 1024;
  time_in[0] = 0;
  time_in[1] = 4;
  time_in[2] = 100;

  err|= setup();
  if(!test_bed->check_for_abort(err)) err|= basic_tests();
  if(my_reader) delete my_reader;
  if(test_contr) delete test_contr;
    
  if(!test_bed->check_for_abort(err)){
    file_prefix = "./files/2dtest/";
    strcpy(block_id, "ey");
    space_in[0] = 0;
    space_in[1] = 1024;
    time_in[0] = 0;
    time_in[1] = 50;
    time_in[2] = 0;

    err|=setup();
    if(!test_bed->check_for_abort(err)) err|= twod_tests();
    if(!test_bed->check_for_abort(err)) err|= twod_space_tests();

  }
  
  test_bed->report_err(err);

  return err;

}

int test_entity_levelone::setup(){
/** \brief Setup to "level one" extraction
*
*
*/

  int err = TEST_PASSED;
  bool use_row_time=false;
  my_reader = new reader(file_prefix, block_id);
  if(my_reader->current_block_is_accum()) use_row_time = true;

  if(!use_row_time){
    n_tims = std::max(time_in[1]-time_in[0], 1);
  }else{
    n_tims = time_in[2];
  }

  size_t n_dims;
  std::vector<size_t> dims;
  int err2 = my_reader->read_dims(n_dims, dims);
  if(err2) err |= TEST_FATAL_ERR;
  
  test_contr = new controller(file_prefix);

  return err;
}

int test_entity_levelone::basic_tests(){
/** \brief Basic tests of process to make levl-1 data
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  int space_dim = space_in[1]-space_in[0];

  dat = data_array(space_dim, n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

  int err2 = my_reader->read_data(dat, time_in, space_in);
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array...
  dat_fft = data_array(space_dim, n_tims);

  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }
  dat.B_ref = get_ref_Bx(file_prefix, space_in, 1);
  //Use first file. 0th accumulator is perhaps broken
  err2 = dat.fft_me(dat_fft);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->set_plasma_B0(dat.B_ref);
  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->generate_spectrum(dat_fft);
  
  int n_dims = dat.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=3){
    lims.push_back(-0.002);
    lims.push_back(0.002);
  }
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);
  
  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
//    dat_fft.write_section_to_file(file, lims);
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}

int test_entity_levelone::twod_tests(){
/** \brief Basic tests of process to make levl-1 data from 2-d input
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  size_t n_dims_in;
  std::vector<size_t> dims_in;
  my_reader->read_dims(n_dims_in, dims_in);
  if(n_dims_in != 2){
    test_bed->report_info("Wrong file dimension", 1);
    return TEST_FATAL_ERR;
  }
  int space_dim = space_in[1]-space_in[0];
//  dat = data_array(space_dim, dims_in[1], n_tims);
  dat = data_array(space_dim, n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

//  int err2 = my_reader->read_data(dat, time_in, space_in);
  int err2 = my_reader->read_data(dat, time_in, space_in, 1);
  //Note this totals on y-dim automagically
  
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array and size FFT accordingly
//  dat_fft = data_array(space_dim, dims_in[1], n_tims);
  dat_fft = data_array();
  dat_fft.clone_empty(dat);
  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }
  dat.B_ref = -1;
  err2 = dat.fft_me(dat_fft);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->make_test_spectrum();
  
  int n_dims = dat.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=3){
    lims.push_back(-0.002);
    lims.push_back(0.002);
  }
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);

  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}
int test_entity_levelone::twod_space_tests(){
/** \brief Basic tests of process to make levl-1 data from 2-d input
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  size_t n_dims_in;
  std::vector<size_t> dims_in;
  my_reader->read_dims(n_dims_in, dims_in);
  if(n_dims_in != 2){
    test_bed->report_info("Wrong file dimension", 1);
    return TEST_FATAL_ERR;
  }
  int space_dim = space_in[1]-space_in[0];
  dat = data_array(space_dim, dims_in[1], n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

  int err2 = my_reader->read_data(dat, time_in, space_in);
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array and size FFT accordingly
//  dat_fft = data_array(space_dim, dims_in[1], n_tims);

  dat.B_ref = -1;

//  dat = dat.total(2);

  dat_fft.clone_empty(dat);
  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }

  err2 = dat.fft_me(dat_fft);

  dat_fft = dat_fft.total(2, 0.01f*my_const.omega_ce, 1.5f*my_const.omega_ce);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->make_test_spectrum();
  
  size_t n_dims = dat_fft.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-0.002);
    lims.push_back(0.002);

  }
  if(n_dims >=3){
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);
  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_k_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}

//----------------------------------------------------------------

test_entity_d::test_entity_d(){

  name = "D checks";
  file_prefix = "./files/";

}
test_entity_d::~test_entity_d(){

  if(test_contr) delete test_contr;
}

int test_entity_d::run(){
/** \todo WRITE!*/
  int err = TEST_PASSED;
  
  deck_constants const_tmp = my_const;
  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
  share_consts();

  test_contr = new controller(file_prefix);
  test_bed->report_info("Reading spectrum", mpi_info.rank);
  //Now dump to file and read back in and compare
  bool err2 = test_contr->add_spectrum("spect_out.dat");
  if(!err2){
    test_bed->report_info("Calculating test D", mpi_info.rank);
    test_contr->add_d(5, 5);
    d_report report = test_contr->get_current_d()->calculate();
    err |= report.error;
    test_bed->report_info("Writing test file", mpi_info.rank);

    std::fstream file;
    file.open("test_d.dat", std::ios::binary|std::ios::trunc|std::ios::out|std::ios::in);
    if(file) test_contr->get_current_d()->write_to_file(file);
    file.close();
  }else{
    err |= TEST_ASSERT_FAIL;
  }
  my_const = const_tmp;
  share_consts();

  test_bed->report_err(err);

  return err;
  
}

//----------------------------------------------------------------

test_entity_bounce::test_entity_bounce(){

}
test_entity_bounce::~test_entity_bounce(){

}

int test_entity_bounce::run(){
/** \todo write!*/
  int err = TEST_PASSED;

  return err;

}

test_entity_nonthermal::test_entity_nonthermal(){
  name = "nonthermal";
}

test_entity_nonthermal::~test_entity_nonthermal(){

}

int test_entity_nonthermal::run(){
/** \brief Check nonthermal module
*
* Check function binding, df/dp values, lookup table
*/
  int err = TEST_PASSED;
  
  //Setup single Max non-thermal
  non_thermal * my_elec = new non_thermal("./files/test");

  calc_type a_par = std::sqrt(2.0)*my_elec->get_v_par() / std::sqrt(1.0 - std::pow(my_elec->get_v_par()/v0, 2));
  calc_type a_perp_sq = 2.0*std::pow(my_elec->get_v_perp(), 2)/(1.0 - 2.0*std::pow(my_elec->get_v_perp()/v0, 2));

  calc_type norm = (my_elec->get_total_dens() - 1.0)/pi/std::sqrt(pi)/a_par/a_perp_sq;

  calc_type df_tmp, f_tmp = my_elec->f_p(0, 0);
  
  if(std::abs(f_tmp / norm -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  const size_t n_trials = 5;
  calc_type inc = 0.1, p_par, p_perp;

  //Test f against expected values
  for(int i=0; i<n_trials*n_trials; i++){
    p_perp = inc* (i % n_trials);
    p_par = inc* (int)(i/n_trials);
    f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
    if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  //Test a couple of extreme cases
  p_perp = 1.5;
  f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
  if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  p_par = 2.5;
  f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
  if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  //Test f against expected values
  for(int i=0; i<n_trials; i++){
    p_perp = inc*i + inc;
    p_par = inc*i + inc;
    //close to 0 deriv -> 0

    f_tmp = my_elec->d_f_p(p_par*a_par, 0, 1);
    if(std::abs(f_tmp /((-2.0)*p_par/a_par*my_elec->f_p(p_par*a_par, 0)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

    f_tmp = my_elec->d_f_p(0, p_perp*std::sqrt(a_perp_sq), 0);
    if(std::abs(f_tmp /((-2.0)*p_perp/std::sqrt(a_perp_sq)*my_elec->f_p(0, p_perp*std::sqrt(a_perp_sq))) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }

  //This next is deep wings so we relax constraint so we're not trying to resolve 1e-60
  df_tmp = my_elec->d_f_p(0, 10.0*std::sqrt(a_perp_sq), 0);
  if(std::abs(df_tmp /((-2.0)*10.0/std::sqrt(a_perp_sq)*my_elec->f_p(0, 10.0*std::sqrt(a_perp_sq))) -1.0) > 2e-6) err |=TEST_WRONG_RESULT;
  
  if(err == TEST_PASSED) my_print("Nonthermal values OK", mpi_info.rank);
  
  err |= test_lookup();
  test_bed->report_err(err);
  return err;

}

int test_entity_nonthermal::test_lookup(){
  /** Check known lookup function against known analytic results. Data generated by generate_lookup in generate_test_data.pro*/
  int err = TEST_PASSED;
  
  non_thermal * my_elec_an = new non_thermal("./files/an");
  non_thermal * my_elec_lookup = new non_thermal("./files/lookup");

  calc_type f_an, f_lookup, p_par, p_perp;
  const size_t n_tests = 100;

  for(size_t i=0; i< n_tests; i++){
    for(size_t j=0; j< n_tests; j++){
      p_par = (float) i *v0/2.0/(float) n_tests;
      p_perp = (float) j *v0/2.0/(float) n_tests;
      f_an = my_elec_an->f_p(p_par, p_perp);
      f_lookup = my_elec_lookup->f_p(p_par, p_perp);
      if(std::abs(f_lookup/f_an - 1.0) > LOW_PRECISION){
        err |= TEST_WRONG_RESULT;
      }
    }
  }
  if(err == TEST_PASSED) my_print("Lookup OK", mpi_info.rank);
  return err;
}

#endif