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
#include "data_array.h"
#include "reader.h"
/*
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
const int TEST_FATAL_ERR = 128;
const int err_tot = 9;
const calc_type PRECISION = 1e-10;/**< Constant for equality at normal precision i.e. from rounding errors etc*/
const calc_type NUM_PRECISION = 1e-6;/**< Constant for equality at good numerical precision, e.g. from numerical integration over 100-1000 pts*/
const calc_type LOW_PRECISION = 5e-3;/**< Constant for equality at low precision, i.e. different approximations to an expression*/
const int max_verbos = 4;
const std::string filename = "tests.log";/**<Test log file*/

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
    std::string get_color_escape(char col);
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
    void set_colour(char col=0);
    bool is_fatal(int err);
    bool check_for_abort(int err);
  
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

bool compare_2d(data_array &lhs, data_array &rhs, bool no_dims_match=false);
bool compare_3d(data_array &lhs, data_array &rhs, bool no_dims_match=false);

/** Test for data array class, assigns values to entry and reads back*/
class test_entity_data_array : public test_entity{
  private:
    data_array test_array;
    int technical_tests();
    int basic_tests();
    int assign();
    int three_d_and_shift();
    int io_tests();
  public:
    test_entity_data_array();
    virtual ~test_entity_data_array();
    virtual int run();

};

/** Combined test: reads test sdf file, stores into data array and runs fft. Test data should be a sine curve with one major frequency which is then checked
*/
class test_entity_get_and_fft : public test_entity{
  private:
    data_array test_dat;
    data_array test_dat_fft;
    reader * test_rdr;
    int one_d();
    int two_d();
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
    int resonant_freq();
    int high_density();
    int other_modes();
    int phi_dom();
    int analytic_dispersion();
  public:
    test_entity_plasma();
    virtual ~test_entity_plasma();
    virtual int run();

};

/** Check spectrum calculations, such as test spectrum derivation etc */
class test_entity_spectrum : public test_entity{
  private:
    data_array test_dat_fft;
    data_array test_spect;
    controller * test_contr;
    std::string file_prefix;
    int setup();
    int basic_tests1();
    int basic_tests2();
    int albertGs_tests();

  public:
    test_entity_spectrum();
    virtual ~test_entity_spectrum();
    virtual int run();
};

/** Full check of deriving a "level one" FFT and spectrum from the various input data formats */
class test_entity_levelone: public test_entity{
  private:
    data_array dat_fft;
    data_array dat;
    controller * test_contr;
    reader * my_reader;
    std::string file_prefix;
    int time_in[3], space_in[2];
    char block_id[ID_SIZE];
    int n_tims;
    int setup();
    int basic_tests();
    int twod_tests();
    int twod_space_tests();
  public:
    test_entity_levelone();
    virtual ~test_entity_levelone();
    virtual int run();
};


/** Spectrum to D test. Setup sample data in a spectrum with analytic solvable form. Calculate resulting D. Cross check*/
class test_entity_d : public test_entity{
  private:
    controller * test_contr;
    std::string file_prefix;

  public:
    test_entity_d();
    virtual ~test_entity_d();
    virtual int run();
};

/** Test bounce averaging. Setup dummy D data across multiple blocks and average. Cross check with analytic results*/
class test_entity_bounce: public test_entity{
  private:

  public:
    test_entity_bounce();
    virtual ~test_entity_bounce();
    virtual int run();
};

class test_entity_nonthermal: public test_entity{
  private:
    int test_lookup();
  public:
    test_entity_nonthermal();
    virtual ~test_entity_nonthermal();
    virtual int run();

};


#endif
#endif