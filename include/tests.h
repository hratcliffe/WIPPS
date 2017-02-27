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

class tests;
extern tests * test_bed; /**< Global testbed, define somewhere in your code*/

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
    bool run_tests();
    void set_verbosity(int verb);
    void set_colour(char col=0);
    bool is_fatal(int err);
    bool check_for_abort(int err);
  
};

#endif
#endif