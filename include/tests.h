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
#include <map>
#include "support.h"
#include <boost/filesystem.hpp>
/*
check 2011 paper for whistler mode tests
Bortnik, 
test particle particle whistler 
quasi-particle treatment of whistlers

*/

/** \defgroup tests Test code
\brief Test framework and tests
*
* Contains the testing framework and information on all tests
*@{
*\defgroup test_bed Test framework
\brief Test framework and constants
*
* Contains test_bed, test entity prototype, constants and testing helper functions
*@{*/

//--- Error code constants to report given sorts of errors ----
const int TEST_PASSED = 0;/**< No error */
const int TEST_WRONG_RESULT = 1;/**< Wrong or unexpected result*/
const int TEST_NULL_RESULT = 2;/**< Null result, no effect or empty object*/
const int TEST_ASSERT_FAIL = 4;/**< Assignment failures, thrown exceptions, invalid results */
const int TEST_USERDEF_ERR1 = 8;/**< Unspecified error, see tests::set_userdef_error()*/
const int TEST_USERDEF_ERR2 = 16;/**< Unspecified error, see tests::set_userdef_error()*/
const int TEST_USERDEF_ERR3 = 32;/**< Unspecified error, see tests::set_userdef_error()*/
const int TEST_USERDEF_ERR4 = 64;/**< Unspecified error, see tests::set_userdef_error()*/
const int TEST_REMOVE_ERR = 128;/**<Temporary flag turning any other non-fatal errors into warnings only*/
const int TEST_FATAL_ERR = 256;/**< Use in conjunction with a specific error if test cannot continue*/
const int err_tot = 10;
const calc_type PRECISION = 1e-10;/**< Constant for equality at normal precision i.e. from rounding errors etc*/
const calc_type NUM_PRECISION = 1e-6;/**< Constant for equality at good numerical precision, e.g. from numerical integration over 100-1000 pts*/
const calc_type LOW_PRECISION = 5e-3;/**< Constant for equality at low precision, i.e. different approximations to an expression*/
const int max_verbos = 4;/**<Verbosity range for messages*/
const std::string filename = "tests.log";/**<Test log file*/
const std::string tests_tmp_dir ="./files/tmp_tests/";/**< Temporary directory for test output.*/


/**\brief Testing instance
*
*Consists of at least a constructor doing any setup required, a name string for output id, a function run() taking no parameters which performs the necessary test and a destructor doing cleanup.
*/
class test_entity{
  public:
    std::string name;/**< The name of the test, which will be reported in the log file*/
    
    test_entity(){;}
    virtual ~test_entity(){;}
    virtual int run()=0;/**<Run the test. Pure virtual means this must be defined in any concrete class @return Error code from test, one of err_codes */

};

/**\brief Test controller
*
*Controls running of tests and their logging etc
*To add a test, do the following:
*
*Descend an object from test_entity which has at least a constructor doing any setup required, a name string for output id, a function run taking no parameters which performs the necessary test and a destructor doing cleanup. Add any other member variables or functions required, including their headers also. In tests::setup_tests create an instance of your class as test_obj = new your_class() and then add your test to the remit using add_test(test_obj); Alternately make the instance and use the global test_bed using test_bed->add(your pntr) from anywhere.
*
*To add errors, add the message into the blank spaces in the list below, err_names, and declare a const int of desired name aliased to TEST_USERDEF_ERR* where * = 1-4
*
*To report the errors by code, call test_bed->report_err(err); To report other salient information use test_bed->report_info(info, verbosity) where the second parameter is an integer describing the verbosity setting at which to print this info (0 = always, larger means more and more detail. 1 is given as default parameter for this function).
*
*To allow conditional running of tests etc, there is a map of runtime_flags which holds all command line arguments. These are assumed to be either a -flag_name or a pair -flag_name integer_flag_value. Available flags are:
*
\htmlinclude tests_runtime_flags.txt
*/
/** To include flags in these docs, include comment line in any functions in the relevant test_entity of the form
"*Set runtime_flag "flag_name" to <action description>*" noting leading and trailing Doxy breaks (* with newline)
*/
/** \todo Remember to include some inverse tests to ensure things are actually doing something \todo Create test data for double-precision!*/

class tests{
  private:
    std::fstream * outfile; /**< Output file handle*/
    int current_test_id;/**< Number in list of test being run*/
    std::vector<test_entity*> test_list;/**< List of tests to run*/
    int verbosity;/**< Verbosity level of output*/
    std::string get_printable_error(int err, int test_id);
    std::string get_colour_escape(char col);
    bool is_fatal(int err);
    void setup_tests();
    void cleanup_tests();
    void create_test_outdir();
    bool no_clean = false;/**< If true, do not clean up testing output directory*/
    void clean_test_outdir(bool no_clean);
  public:
    std::map<std::string, int> runtime_flags;/**<Strings for runtime control*/
    tests(){setup_tests();};
    ~tests(){cleanup_tests();};
    void report_err(int err, int test_id=-1);
    void report_info(std::string info, int verb_to_print = 1);
    void add_test(test_entity* test);
    bool run_tests();
    void set_verbosity(size_t verb);
    void set_runtime_flags(int argc, char *argv[]);
    void set_colour(char col=0);
    bool check_for_abort(int err);
    bool set_userdef_error(int err_code, std::string message);
};

extern tests * test_bed; /**< Global testbed, defined in tests.cpp*/

/** @} */
/** @} */

#endif
#endif