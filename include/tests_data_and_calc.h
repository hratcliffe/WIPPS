//
//  tests_data_and_calc.h
//  
//
//  Breakout by Heather Ratcliffe on 27/02/2017.
//
//
#ifdef RUN_TESTS_AND_EXIT

#ifndef _tests_data_h
#define _tests_data_h

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "tests.h"
#include "support.h"
#include "data_array.h"
#include "reader.h"
#include "controller.h"

class plasma;

/** \ingroup tests
*\defgroup tests_data Calculation tests
\brief Test physics and calculations
*
* Tests for physics modelling and calculation details
*@{*/

my_type calc_I_omega(my_type omega, spectrum * my_spect, controller * my_contr);

/**Check plasma functions, get_omega and dispersion relation
*/
class test_entity_plasma : public test_entity{
  private:
    plasma * plas;/**< Pointer to plasma object under test*/
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
    data_array test_dat_fft;/**< Holds FFTd test data*/
    data_array test_spect;/**< Holds test spectral data*/
    controller * test_contr;/**< Controller for testing*/
    std::string file_prefix;/**< File prefix prepended to all test files*/
    int setup();
    int basic_tests1();
    int basic_tests2();
    int albertGs_tests();
    int technical_tests();

  public:
    test_entity_spectrum();
    virtual ~test_entity_spectrum();
    virtual int run();
};

/** Full check of deriving a "level one" FFT and spectrum from the various input data formats */
class test_entity_levelone: public test_entity{
  private:
    data_array dat_fft;/**< Holds FFTd test data*/
    data_array dat;/**< Holds raw test data*/
    controller * test_contr;/**< Controller for testing*/
    reader * my_reader;/**< File reader for testing*/
    std::string file_prefix;/**< File prefix prepended to all test files*/
    size_t time_in[3];/**< Time specs for data to read*/
    size_t space_in[2];/**< Space specs for data to read*/
    char block_id[ID_SIZE];/**< Name of block to read*/
    int n_tims;/**<Number of times to read*/
    int setup();
    int basic_tests(size_t n_dims_in, int flatten_on, bool has_freq, std::string outfile_tag = "", int total_fft = -1, my_type band_min = 0.0, my_type band_max = 0.0);
  public:
    test_entity_levelone();
    virtual ~test_entity_levelone();
    virtual int run();
};


/** Spectrum to D test. Setup sample data in a spectrum with analytic solvable form. Calculate resulting D. Cross check*/
class test_entity_d : public test_entity{
  private:
    controller * test_contr;/**< Controller for testing*/
    std::string file_prefix;/**< File prefix prepended to test files*/
    int basic_tests();
    int full_D_tests();
    data_array read_padie_data(bool single_n, int n);
  public:
    test_entity_d();
    virtual ~test_entity_d();
    virtual int run();
};

/** Test bounce averaging. Setup dummy D data across multiple blocks and average. Cross check with analytic results*/
class test_entity_bounce: public test_entity{
  private:
    controller * test_contr;/**< Controller for testing*/
    std::string file_prefix;/**< File prefix prepended to test files*/
    int bounce_cases(bounce_av_data bounce_dat);
  public:
    test_entity_bounce();
    virtual ~test_entity_bounce();
    virtual int run();
};

/** Test non-thermal electron specification*/
class test_entity_nonthermal: public test_entity{
  private:
    int test_lookup();
  public:
    test_entity_nonthermal();
    virtual ~test_entity_nonthermal(){;}
    virtual int run();

};

/** @} */

#endif
#endif