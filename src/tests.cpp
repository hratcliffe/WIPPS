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
#include <unistd.h>
#include <functional>
#ifdef USE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
#endif
#include "tests.h"
#include "tests_basic_and_code.h"
#include "tests_data_and_calc.h"

tests * test_bed;

/** \ingroup test_bed */
const int err_codes[err_tot] ={TEST_PASSED, TEST_WRONG_RESULT, TEST_NULL_RESULT, TEST_ASSERT_FAIL, TEST_USERDEF_ERR1, TEST_USERDEF_ERR2, TEST_USERDEF_ERR3, TEST_USERDEF_ERR4, TEST_REMOVE_ERR, TEST_FATAL_ERR};/**< List of error codes available*/

/** \ingroup test_bed */
std::string err_names[err_tot]={"None", "Wrong result", "Invalid Null result", "Assignment or assertion failed", "{Message 1 here}", "{Message 2 here}", "{Message 3 here}", "{Message 4 here}", "SQUASHED ERRORS", "Fatal error"};/**< Names corresponding to error codes, which are reported in log files*/

void tests::create_test_outdir(){
/** \brief Create the test directory
*
* Creates the directory tests_tmp_dir for temporary test output. If the directory exists, nothing is done and the no_clean flag is set
*/
#ifdef USE_BOOST_FILESYSTEM
  boost::filesystem::path dir(tests_tmp_dir.c_str());
  if(boost::filesystem::exists(dir)){
    no_clean = true;
  }else{
    this->report_info("Creating temporary dir "+tests_tmp_dir, 0);
    boost::filesystem::create_directory(dir);
  }
#endif
}

void tests::clean_test_outdir(bool no_clean){
/** \brief Clean the test directory
*
* Removes contents and directory tests_tmp_dir, unless the no_clean parameter is set, when nothing is done. Note there may be other flags to consider so this is not as stupid as it seems. NB this writes to log file!
@param no_clean If set, do nothing
*/
#ifdef USE_BOOST_FILESYSTEM
  if(!no_clean){
    this->report_info("Removing temporary dir "+tests_tmp_dir, 0);
    boost::filesystem::path dir(tests_tmp_dir.c_str());
    boost::filesystem::remove_all(dir);
  }
#endif
}

void tests::set_verbosity(size_t verb){
/** \brief Set verbosity
*
*Set the verbosity of testing output, from 0 (minimal) to max_verbos cap
@param verb Verbosity level to set
*/
  this->verbosity = std::min((int)verb, max_verbos);
}

void tests::set_runtime_flags(int argc, char *argv[]){
/** Set the runtime flags
*
* Takes the command line args and stores in the runtime flags map. Each flag can be either a -thing or a -thing val where val is an int. Flags may not start with a digit
@param argc Command line parameter number
@param argv Command line parameter list
*/

  std::string current_flag;
  int val = 0;
  int i = 1;
  while(i < argc){
    current_flag = argv[i];
    if(argv[i][0] == '-') current_flag = current_flag.substr(1, current_flag.size());
    else continue;
    val = 0;
    if(i < argc-1 && (argv[i+1][0] != '-'  || (argv[i+1][1] >='0' && argv[i+1][1] <='9'))){
      val = checked_strtol(argv[i+1]);
      i++;
    }
    runtime_flags[str_to_lower(current_flag)] = val;
    i++;
  }
}

void tests::add_test(test_entity * test){
/** Adds a test to the list to be performed
@param test Pointer to test object to add. Object will be deleted on completion */
  test_list.push_back(test);
  report_info("  Added "+test->name, 2);
}

bool tests::run_tests(){
/** \brief Run scheduled tests
*
*Runs each test in list and reports total errors found. 
@return 0 for no errors, 1 else
*/

  int total_errs = 0, total_warnings = 0;
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
    int err = test_list[current_test_id]->run();
    report_err(err);
    if(err != TEST_PASSED && err != TEST_REMOVE_ERR){
      if((err & TEST_REMOVE_ERR) != TEST_REMOVE_ERR) total_errs++;
      else total_warnings ++;
    }
    //Add one if any error was returned
  }
  this->set_colour('*');
  if(total_errs > 0){
    this->set_colour('r');
  }else if(total_warnings > 0){
    this->set_colour('y');
  }else{
    this->set_colour('b');
  }
  my_error_print(mk_str(total_errs)+" failed tests", mpi_info.rank);
  
  if(total_warnings > 0){
    this->set_colour('y');
    my_error_print(mk_str(total_warnings)+" tests with squashed errors", mpi_info.rank);
  }
  this->set_colour();
  return total_errs > 0;

}

bool tests::is_fatal(int err){
/** \brief Checks whether error is fatal
*
* @param err Error code, one of err_codes
@return Boolean true if fatal, false else
*/

  if((err & TEST_FATAL_ERR) == TEST_FATAL_ERR) return 1;
  else return 0;
}

bool tests::check_for_abort(int err){
/**\brief Check for abort condition
*
*Checks whether current test has had fatal error
@param err Error param to check
@return true if test should abort, false else
*/
  if(is_fatal(err)){
    set_colour('r');
    set_colour('*');
    my_error_print("Fatal error occured. Aborting test "+test_list[current_test_id]->name, mpi_info.rank);
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
* Logs error text corresponding to code err for test defined by test_id.
@param err Error to report
@param test_id Index of corresponding test
*/
  if(test_id == -1) test_id = current_test_id;
  if(err ==TEST_PASSED || err == TEST_REMOVE_ERR){
    set_colour('b');
  }else if((err & TEST_REMOVE_ERR) == TEST_REMOVE_ERR){
    set_colour('y');
  }else{
    set_colour('r');
  }
  if(is_fatal(err)) set_colour('*');
  my_error_print(outfile, get_printable_error(err, test_id), mpi_info.rank);
  my_error_print(nullptr, get_printable_error(err, test_id), mpi_info.rank);
  set_colour();

}

void tests::report_info(std::string info, int verb_to_print){
/** \brief Other test info
*
*Records string info to the tests.log file and to screen, according to requested verbosity.
@param info Information to print
@param verb_to_print Verbosity level above which to print this info
*/
  if(verb_to_print <= this->verbosity){
    my_print(outfile, info, mpi_info.rank);
    my_print(nullptr, info, mpi_info.rank);
  
  }

}

std::string tests::get_printable_error(int err, int test_id){
/** \brief Make an error message
*
* Converts error code to printable string, adds code for reference and adds test name. Note code is bitmask and additional errors are appended together
@param err Error code
@param test_id Index of current test in list
@return Error information as printable string
*/
  std::string err_string="";
  int err_remaining = err;
  if(err_remaining != TEST_PASSED){
    for(int i=err_tot-1; i>0; --i){
      //Run most to least significant
      if((err_remaining & err_codes[i]) == err_codes[i]){
        err_string +=err_names[i];
        err_remaining -= err_codes[i];
        if(err_remaining != TEST_PASSED) err_string +=" &";
        err_string += " ";
      }
    }
  
    err_string = "Error "+err_string+"(code "+mk_str(err)+") on";
  }
  else err_string = "Passed";
  return err_string+" test "+test_list[test_id]->name;

}

/**
     * @class dummy_colour
     * Available colour codes are RGB, CMYK and W(hite) plus * (bold) _ (underlined) ? (blink, irritating and not supported on some terminals) and $ (reverse foreground and background). A call with no arguments, or with 0 or '0' reset to default.
*/

void tests::set_colour(char col){
/** \brief Set output text colour
*
*Set terminal output colour using std escape sequences. NB technically not MPI safe. Use sparingly to highlight important information. \copydoc dummy_colour
@param col Colour code, one of rgb cmyk or [w]hite, or format code,*_?$ 0, '0' or blank
*/
  if(isatty(fileno(stdout))) my_print(this->get_colour_escape(col), mpi_info.rank, 0, true);
}

inline std::string tests::get_colour_escape(char col){
/** \brief
*This returns the terminal escape string to set given colour.
*
\copydoc dummy_colour
@param col The colour code
@return String giving VT100 escape for colour
*/
  if(col >='A' and col <='Z') col += 32;
  //ASCII upper to lower
  switch (col) {
    case 0:
    case '0':
      return "\033[0m";
      break;//"break" is Redundant but clearer
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

/* Add tests here!*/
void tests::setup_tests(){
/** \brief Setup test bed
*
*Opens reporting file, creates temporary dir. Then instantiates all the test objects and adds them into the test_list
\todo Clean up error codes
\todo Clean up printing verbosities
*/
  set_verbosity(max_verbos);//"All" info is printed during setup

  outfile = new std::fstream();
  outfile->open(filename.c_str(), std::ios::out);
  if(!outfile->is_open()){
    my_error_print("Error opening "+filename, mpi_info.rank);
    //can't log so return with empty test list
    return;
  }

  create_test_outdir();

  test_entity * test_obj;

  //Assign first user-def error to be "File IO fail"
  set_userdef_error(TEST_USERDEF_ERR1, "File IO failed");

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
  test_obj = new test_entity_levelone();
  add_test(test_obj);
  test_obj = new test_entity_d();
  add_test(test_obj);
  test_obj = new test_entity_nonthermal();
  add_test(test_obj);
  test_obj = new test_entity_bounce();
  add_test(test_obj);

}

void tests::cleanup_tests(){
/** \brief Clean up after testing
*
*Deletes test objects, and closes logfile
*/
  if(runtime_flags.count("no_clean") == 0){
    this->clean_test_outdir(this->no_clean);
  }

  if(outfile->is_open()){
    this->report_info("Testing complete and logged in " +filename, 0);
    outfile->close();
  }else{
    this->report_info("No logfile generated", 0);

  }
  delete outfile;
  
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
    delete test_list[current_test_id];
    test_list[current_test_id] = nullptr;
  }
}

bool tests::set_userdef_error(int err_code, std::string message){
/** \brief Set a user-defined error string
*
*Sets message associated with err_code to message string. NOTE this affects code associated with past and future use of the err_code string. 
@param err_code Code to set, one of the TEST_USERDEF_ERR[x]
@param message Message for this error
@return 0 success, 1 for invalid change (trying to change a non-user error)
*/

  //We assume ordering of the USERDEF errors, but nothing more
  //Check for out-of-range
  if(err_code < TEST_USERDEF_ERR1 || err_code > TEST_USERDEF_ERR4) return 1;
  //Find index and update message
  for(size_t i = 0; i<= err_tot; i++){
    if(err_codes[i] == err_code){
      err_names[i] = message;
      return 0;
    }
  }
  return 1;//Was invalid code
}

#endif

