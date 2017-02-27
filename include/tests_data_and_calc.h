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

class plasma;
class controller;

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
    int technical_tests();

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
    size_t time_in[3];
    size_t space_in[2];
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

/** Test non-thermal electron specification*/
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