//
//  tests.h
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//
#ifdef RUN_TESTS_AND_EXIT

#ifndef _tests_h
#define _tests_h

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "support.h"

/*
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
const calc_type PRECISION = 1e-10;/**< Constant for equality at normal precision i.e. from rounding errors etc*/
const calc_type LOW_PRECISION = 5e-3;/**< Constant for equality at low precision, i.e. different approximations to an expression*/
const int max_verbos = 4;
const std::string filename = "tests.log";/**<Test log file*/

class reader;
class data_array;
class plasma;
class controller;

/**\brief Testing instance
*
*Consists of at least a constructor doing any setup required, a name string for output id, a function run() taking no parameters which performs the necessary test and a destructor doing cleanup.
*/
class test_entity{
public:
  std::string name;/**< The name of the test, which will be reported in the log file*/
  
  test_entity(){;}
  virtual ~test_entity(){;}
  virtual int run()=0;/*Pure virtual because we don't want an instances of this template*/

};

/**\brief Test controller
*
*Controls running of tests and their logging etc
*To add a test, do the following:
*Descend an object from test_entity which has at least a constructor doing any setup required, a name string for output id, a function run taking no parameters which performs the necessary test and a destructor doing cleanup. Add any other member variables or functions required, including their headers also. In tests::setup_tests create an instance of your class as test_obj = new your_class() and then add your test to the remit using add_test(test_obj); Alternately make the instance and use the global test_bed using test_bed->add(your pntr) from anywhere.
*To add errors, add the message into the blank spaces in the list below, err_names, and declare a const int of desired name aliased to TEST_USERDEF_ERR* where * = 1-4
*To report the errors by code, call test_bed->report_err(err); To report other salient information use test_bed->report_info(info, verbosity) where the second parameter is an integer describing the verbosity setting at which to print this info (0=always, the larger int means more and more detail).


*/

class tests{
private:

  std::string get_printable_error(int err, int test_id);
  std::fstream * outfile; /**< Output file handle*/
  int current_test_id;/**< Number in list of test being run*/
  std::vector<test_entity*> test_list;/**< List of tests to run*/
  int verbosity;/**< Verbosity level of output*/
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

/** Test for reader class */
class test_entity_reader : public test_entity{
  private:
  reader * test_rdr;
  reader * accum_reader;
  const static int size = 49367784;
  //Size of my test file...

  public:
  test_entity_reader();
  virtual ~test_entity_reader();
  virtual int run();

};

/** Test for data array class, assigns values to entry and reads back*/
class test_entity_data_array : public test_entity{
  private:
  data_array * test_array;
  public:
  test_entity_data_array();
  virtual ~test_entity_data_array();
  virtual int run();

};

/** Combined test: reads test sdf file, stores into data array and runs fft. Test data should be a sine curive with one major frequency which is then checked
*/
class test_entity_get_and_fft : public test_entity{
  private:
  data_array * test_dat;
  data_array * test_dat_fft;
  reader * test_rdr;
  public:
  test_entity_get_and_fft();
  virtual ~test_entity_get_and_fft();
  virtual int run();

};

/** Test basic maths routines, including where, integrator, boxcar smoothing, interpolation and cubic solver. Runs predefined test problems and tests results.
*/
class test_entity_basic_maths : public test_entity{
  private:
    calc_type * data_square;
    calc_type * data_positive;
    calc_type * data_tmp;
    calc_type * axis;
    my_type * axisf;
    calc_type * d_axis;
    int size;
  public:
  test_entity_basic_maths();
  virtual ~test_entity_basic_maths();
  virtual int run();

};

/** Test external maths routines use and interpretation etc. Currently Bessel functions from boost
*/
class test_entity_extern_maths : public test_entity{
  private:

  public:
  test_entity_extern_maths();
  virtual ~test_entity_extern_maths();
  virtual int run();

};

/**Check plasma functions, get_omega and dispersion relation
*/
class test_entity_plasma : public test_entity{
  private:
  plasma * plas;
  public:
  test_entity_plasma();
  virtual ~test_entity_plasma();
  virtual int run();

};

/** Check basic spectrum calculations, such as test spectrum derivation etc */
class test_entity_spectrum : public test_entity{
  private:
  data_array * test_dat_fft;
  data_array * test_spect;
  controller * test_contr;
  std::string file_prefix;
  int tim_in[3], space_in[2];

  public:
  test_entity_spectrum();
  virtual ~test_entity_spectrum();
  virtual int run();
  int setup();
  int basic_tests();
  int albertGs_tests();
};

/**Check G1 ad G2 from Albert \todo Write
*/

#endif
#endif