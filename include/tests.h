//
//  tests.h
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

#ifndef _tests_h
#define _tests_h

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "support.h"

//OK this will contain all our basic code testing, interfacing with libraries etc etc
//Maths tests can be seperate
/** 
[10/11/2015 15:24:26] Heather Ratcliffe: Hmm, OK so for my code what I mainly want is to have a test .sdf file with a couple of blocks in, and to read from it say a 1-d and 2-d array and check they match what I put in them
[10/11/2015 15:25:01] Christopher Brady: Yep. Probably just generating the file from EPOCH or LARE is the easiest solution.
[10/11/2015 15:25:22] Heather Ratcliffe: Similarly generate 1-d and 2-d arrays with sine waves, and axes, and check the FFTs give the right result and frequency
[10/11/2015 15:25:52] Christopher Brady: Yep. Could combine the two and just adjust the arrays in EPOCH or LARE just before writing the file.
[10/11/2015 15:27:10] Heather Ratcliffe: Basic tests of numerical integrals and the like
[10/11/2015 15:27:52] Heather Ratcliffe: And some evil evil tests of the dispersion solvers and the like and of the overall coefficient production

check 2011 paper for whistler mode tests
Bortnik, 
test particle particle whistler 
quasi-particle treatment of whistlers

*/

const int TEST_PASSED = 0;
const int TEST_WRONG_RESULT = 1;
const int TEST_NULL_RESULT = 2;
const int TEST_ASSERT_FAIL = 4;
const int TEST_USERDEF_ERR1 = 8;
const int TEST_USERDEF_ERR2 = 16;
const int TEST_USERDEF_ERR3 = 32;
const int TEST_USERDEF_ERR4 = 64;
const int err_tot = 8;
const calc_type PRECISION = 1e-6;
const int max_verbos = 1;
const std::string filename = "tests.log";

class reader;
class data_array;

class test_entity{
//Consists of at least a constructor doing any setup required, a name string for output id, a function run taking no parameters which performs the necessary test and a destructor doing cleanup.
public:
  std::string name;
  
  test_entity(){;}
  virtual ~test_entity(){;}
  virtual int run()=0;
  //Pure virtual because we don't want an instances of this template

};

class tests{
private:
  int test_template();

  std::string get_printable_error(int err, int test_id);
  std::fstream * outfile;
  int current_test_id;
  std::vector<test_entity*> test_list;
  int verbosity;
public:
  void report_err(int err, int test_id=-1);
  void report_info(std::string info, int verb_to_print = 1, int test_id=-1);
  tests();
  ~tests();
  void setup_tests();
  void add_test(test_entity* test);
  void cleanup_tests();
  void run_tests();
  void set_verbosity(int verb);
  
};

class test_entity_reader : public test_entity{
  private:
  reader * test_rdr;
  const static int size = 49367784;
  //Size of my test file...

  public:
  test_entity_reader();
  virtual ~test_entity_reader();
  virtual int run();

};
class test_entity_data_array : public test_entity{
  private:
  data_array * test_array;
  public:
  test_entity_data_array();
  virtual ~test_entity_data_array();
  virtual int run();

};
class test_entity_get_and_fft : public test_entity{
  private:
  data_array * test_dat;
  data_array * test_dat_fft;
  reader * test_rdr;
  public:
  test_entity_get_and_fft();
  virtual ~test_entity_get_and_fft();
  virtual int run();

//A chunky test. It should read some test data which will be a sine curve from an "ex" block in a test sdf file, into a data array, and then fft it. We then check that the resulting major frequency is as expected
};
class test_entity_basic_maths : public test_entity{
  private:
    calc_type * data_square;
    calc_type * data_positive;
    calc_type * data_tmp;
    calc_type * axis;
    calc_type * d_axis;
    int size;
  public:
  test_entity_basic_maths();
  virtual ~test_entity_basic_maths();
  virtual int run();

//Verify the basic maths functions after any changes. E.g. integrator, boxcar smooth etc
};

class test_entity_extern_maths : public test_entity{
  private:

  public:
  test_entity_extern_maths();
  virtual ~test_entity_extern_maths();
  virtual int run();

//check use and interpretation of external maths, such as boost Bessel fns etc
};


#endif
