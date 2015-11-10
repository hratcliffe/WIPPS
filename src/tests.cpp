//
//  tests.cpp
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

#include <stdio.h>
#include "tests.h"
#include "reader.h"

class reader;
void tests::setup_tests(){
/** \brief Generate necessary objects for tests
*
*
*/

  char block_id[10]="ex";
  test_reader = new reader("test", block_id);


}

void tests::cleanup_tests(){
/** \brief Delete test objects
*
*
*/

  delete test_reader;

}

void tests::run_tests(){
/** \brief Delete test objects
*
*
*/
  int err = TEST_PASSED;

  err = test_reader();


}

int tests::test_reader(){
  int err = TEST_PASSED;

/* do testing */

  report_err(err);
  return err;

}
