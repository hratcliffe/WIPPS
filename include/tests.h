//
//  tests.h
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//

#ifndef _tests_h
#define _tests_h

//OK this will contain all our basic code testing, interfacing with libraries etc etc
//Maths tests can be seperate
/** 
[10/11/2015 15:24:26] Heather Ratcliffe: Hmm, OK so for my code what I mainly want is to have a test .sdf file with a couple of blocks in, and to read from it say a 1-d and 2-d array and check they match what I put in them
[10/11/2015 15:25:01] Christopher Brady: Yep. Probably just generating the file from EPOCH or LARE is the easiest solution.
[10/11/2015 15:25:22] Heather Ratcliffe: Similarly generate 1-d and 2-d arrays with sine waves, and axes, and check the FFTs give the right result and frequency
[10/11/2015 15:25:52] Christopher Brady: Yep. Could combine the two and just adjust the arrays in EPOCH or LARE just before writing the file.
[10/11/2015 15:26:07] Heather Ratcliffe: Yeah I thought about that
[10/11/2015 15:26:55] Heather Ratcliffe: Couple of quick more traditional tests that the data_array class is storing the data properly, and it’s methods to extract and populate elements work
[10/11/2015 15:26:59] Christopher Brady: It's a bit tricky to write SDF files otherwise.
[10/11/2015 15:27:08] Christopher Brady: Yep. Makes sense.
[10/11/2015 15:27:10] Heather Ratcliffe: Basic tests of numerical integrals and the like
[10/11/2015 15:27:33] Christopher Brady: Yep, good approach.
[10/11/2015 15:27:52] Heather Ratcliffe: And some evil evil tests of the dispersion solvers and the like and of the overall coefficient production
*/

class reader;

const int TEST_PASSED = 0;
const int TEST_WRONG_RESULT = 1;
const int TEST_NULL_RESULT = 2;


class tests{
private:
  int test_reader();

public:
  reader * test_reader;

void setup_tests();
void cleanup_tests();
void run_tests();










};



#endif
