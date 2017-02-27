//
//  tests_basic_and_code.h
//  
//
//  Breakout by Heather Ratcliffe on 27/02/2017.
//
//
#ifdef RUN_TESTS_AND_EXIT

#ifndef _tests_basic_h
#define _tests_basic_h

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "tests.h"
#include "support.h"
#include "data_array.h"
#include "reader.h"


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

bool compare_2d(data_array const &lhs, data_array const &rhs, bool no_dims_match=false);
bool compare_3d(data_array &lhs, data_array &rhs, bool no_dims_match=false);

/** Test for data array class, assigns values to entry and reads back*/
class test_entity_data_array : public test_entity{
  private:
    data_array test_array;
    my_type total_all;
    int technical_tests();
    int basic_tests();
    int assign();
    int three_d_and_shift();
    int io_tests();
    int set_vals_2d();
    int set_vals_3d();
  public:
    test_entity_data_array();
    virtual ~test_entity_data_array(){;}
    virtual int run();

};

/** Combined test: reads test sdf file, stores into data array and runs fft. Test data should be a sine curve with one major frequency which is then checked
*/
class test_entity_get_and_fft : public test_entity{
  private:
    reader * test_rdr;
    int one_d();
    int two_d();
    int fft_and_check_1d(data_array & dat_in, data_array & dat_fft, my_type expected_max, bool single_max = false);
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
    const int size = 256;//Size of arrays to use for maths checks
    void setup_arrays();
    void teardown_arrays();
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


#endif
#endif