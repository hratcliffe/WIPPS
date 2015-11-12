//
//  tests.cpp
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

/** \file tests.cpp 
*
*To add a test, do the following:
*Descend an object from test_entity which has at least a constructor doing any setup required, a name string for output id, a function run taking no parameters which performs the necessary test and a destructor doing cleanup. Add any other member variables or functions required. In tests::setup_tests create an instance of your class as test_obj = new your_class() and then add your test to the remit using add_test(test_obj); Alternately make the instance and use the global test_bed using test_bed->add(your pntr) from anywhere

*/

#include <stdio.h>
#include "tests.h"
#include "reader.h"
#include "support.h"

extern mpi_info_struc mpi_info;
extern tests * test_bed;
class reader;


const int TEST_PASSED = 0;
const int TEST_WRONG_RESULT = 1;
const int TEST_NULL_RESULT = 2;
const std::string filename = "tests.log";

std::string err_names[3]={"None", "Wrong result", "Invalid Null result"};

tests::tests(){
  setup_tests();
}
tests::~tests(){
  cleanup_tests();
}

void tests::setup_tests(){
/** \brief Setup test bed
*
*Opens reporting file. Then instantiates all the test objects and adds them into the test_list
*/
  outfile = new std::fstream();
  outfile->open(filename.c_str(), std::ios::out);
  if(!outfile->is_open()){
    std::cout<<"Error opening "<<filename<<std::endl;
    //can't log so return with empty test list
    return;
  }

  test_entity * test_obj;

  test_obj = new test_entity_reader();
  add_test(test_obj);
  //these two lines are needed for each test you want to do.

}
void tests::add_test(test_entity * test){
  /** Adds a test to the list to be performed*/
  test_list.push_back(test);
}

void tests::report_err(int err, int test_id){
  if(test_id == -1) test_id = current_test_id;

  my_print(outfile, get_printable_error(err, test_id), mpi_info.rank);

}

std::string tests::get_printable_error(int err, int test_id){

  std::string err_string;
  if(err!=TEST_PASSED) err_string = "Error "+err_names[err]+" on";
  else err_string = "Passed";
  return err_string+" test "+test_list[test_id]->name;

}

void tests::cleanup_tests(){
/** \brief Delete test objects
*
*
*/
  if(outfile->is_open()){
    std::cout<<"Testing complete and logged in "<<filename<<std::endl;
    outfile->close();
  }else{
    std::cout<<"No log generated"<<std::endl;
  }
  delete outfile;
  for(current_test_id=0; current_test_id<test_list.size(); current_test_id++){
    delete test_list[current_test_id];
  
  }

  
}

void tests::run_tests(){
/** \brief Delete test objects
*
*
*/
  int err = TEST_PASSED;
  for(current_test_id=0; current_test_id<test_list.size(); current_test_id++){
    test_list[current_test_id]->run();
    
  }

}

test_entity_reader::test_entity_reader(){
  name = "reader class";
  test_rdr = new reader("test", "run_info");
  
}
test_entity_reader::~test_entity_reader(){

  delete test_rdr;
}

int test_entity_reader::run(){
  int err = TEST_PASSED;
  
  int sz = test_rdr->get_file_size();
  if(sz != size) err |= TEST_WRONG_RESULT;

/* do testing */
  test_bed->report_err(err);
  return err;

}

test_entity_data_array::test_entity_data_array(){
  name = "data array class";
  
}
test_entity_data_array::~test_entity_data_array(){

}

int test_entity_data_array::run(){
  int err = TEST_PASSED;
  

/* do testing */
  test_bed->report_err(err);
  return err;

}

