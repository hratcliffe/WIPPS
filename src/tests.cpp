//
//  tests.cpp
//  
//
//  Created by Heather Ratcliffe on 10/11/2015.
//
//


#include <stdio.h>
#include <math.h>
#include <cmath>
#include "tests.h"
#include "reader.h"
#include "support.h"
#include "plasma.h"
#include "my_array.h"
#include <math.h>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more


extern tests * test_bed; /**< Global testbed, define somewhere in your code*/

extern mpi_info_struc mpi_info;
extern deck_constants my_const;

const int err_codes[err_tot] ={TEST_PASSED, TEST_WRONG_RESULT, TEST_NULL_RESULT, TEST_ASSERT_FAIL, TEST_USERDEF_ERR1, TEST_USERDEF_ERR2, TEST_USERDEF_ERR3, TEST_USERDEF_ERR4};/**< List of error codes available*/

std::string err_names[err_tot]={"None", "Wrong result", "Invalid Null result", "Assignment or assertion failed", "", "", "", ""};/**< Names corresponding to error codes, which are reported in log files*/

tests::tests(){
  setup_tests();
}
tests::~tests(){
  cleanup_tests();
}

void tests::set_verbosity(int verb){
/** Set the verbosity of testing output, from 0 (minimal) to max_verbos.*/
  if((verb > 0)) this->verbosity = std::max(verb, max_verbos);

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

}
void tests::add_test(test_entity * test){
  /** Adds a test to the list to be performed*/
  test_list.push_back(test);
}

void tests::report_err(int err, int test_id){
/** \brief Log error
*
* Logs error text corresponding to code err for test defined by test_id. Errors are always recorded.*/
  if(test_id == -1) test_id = current_test_id;

  my_print(outfile, get_printable_error(err, test_id), mpi_info.rank);
  my_print(nullptr, get_printable_error(err, test_id), mpi_info.rank);

}

void tests::report_info(std::string info, int verb_to_print, int test_id){
/** \brief Other test info
*
*Records string info to the tests.log file and to screen, according to requested verbosity.
*/
  if(test_id == -1) test_id = current_test_id;
  if(verb_to_print <= this->verbosity){
    my_print(outfile, info, mpi_info.rank);
    my_print(nullptr, info, mpi_info.rank);
  
  }

}

std::string tests::get_printable_error(int err, int test_id){
/** \brief Make an error message
*
* Converts error code to printable string, adds code for reference and adds test name. Note code is bitmask and additional errors are appended together
*/
  std::string err_string="";
  if(err!=TEST_PASSED){
    for(int i=err_tot-1; i>0; --i){
      //Run most to least significant
      if((err & err_codes[i]) == err_codes[i]){
        err_string +=err_names[i] + ", ";
      }
    }
  
    err_string = "Error "+err_string+"(code "+mk_str(err)+") on";
  }
  else err_string = "Passed";
  return err_string+" test "+test_list[test_id]->name;

}

void tests::cleanup_tests(){
/** \brief Delete test objects
*
*
*/
  if(outfile->is_open()){
    this->report_info("Testing complete and logged in " +filename, 0);
    outfile->close();
  }else{
    this->report_info("No logfile generated", 0);

  }
  delete outfile;
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
    delete test_list[current_test_id];
  
  }

  
}

void tests::run_tests(){
/** \brief Run scheduled tests
*
*
*/
  for(current_test_id=0; current_test_id< (int)test_list.size(); current_test_id++){
    test_list[current_test_id]->run();
    
  }

}

test_entity_reader::test_entity_reader(){
  name = "reader class";
  char block_id[10]= "run_info";
  test_rdr = new reader("test", block_id);
  
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
  test_array = new data_array(10, 10);
  
}
test_entity_data_array::~test_entity_data_array(){

  delete test_array;
}

int test_entity_data_array::run(){
  int err = TEST_PASSED;
  bool tmp_err;
  
  int val;
  if(test_array->get_dims() == 2){
    //assign each element to unique val

    for(int i=0; i<test_array->get_dims(0); i++){
      for(int j =0; j<test_array->get_dims(1); j++){
        tmp_err=test_array->set_element(i, j, (i+1)*(2*j+1));
        if(tmp_err) err |= TEST_ASSERT_FAIL;
      }
    }

    //test assignments worked

    for(int i=0; i<test_array->get_dims(0); i++){
      for(int j =0; j<test_array->get_dims(1); j++){
        val = test_array->get_element(i,j);
        if(val != (i+1)*(2*j+1)) err |=TEST_WRONG_RESULT;
      }
    }
  }

/* do testing */
  test_bed->report_err(err);
  return err;

}

test_entity_get_and_fft::test_entity_get_and_fft(){
  name = "read and FFT";
  char block_id[10]= "ex";
  test_rdr = new reader("test", block_id);

}
test_entity_get_and_fft::~test_entity_get_and_fft(){
  delete test_rdr;
  delete test_dat;
  delete test_dat_fft;
  

}

int test_entity_get_and_fft::run(){

  int err = TEST_PASSED;

  int tim_in[2], space_in[2];
  tim_in[0]=0;
  tim_in[1]=1;
  space_in[0]=0;
  space_in[1]=-1;

  int n_tims = std::max(tim_in[1]-tim_in[0], 1);

  int n_dims;
  std::vector<int> dims;
  test_rdr->read_dims(n_dims, dims);
//  std::cout<<dims[0]<<std::endl;
  if(n_dims !=1){
    err |= TEST_WRONG_RESULT;
    if(err == TEST_PASSED) test_bed->report_info("Array dims wrong", 1);
    test_bed->report_err(err);

    return err;
    //nothing more worth doing right now...
  }

  test_dat = new data_array(dims[0], n_tims);
  test_dat_fft = new data_array(dims[0], n_tims);
  if(!test_dat->is_good()||!test_dat_fft->is_good()){
    err|=TEST_ASSERT_FAIL;
    return err;
  }

  test_rdr->read_data(test_dat, tim_in, space_in);

  bool tmp_err = test_dat->fft_me(test_dat_fft);
  if(tmp_err) err|=TEST_ASSERT_FAIL;
  if(err == TEST_PASSED) test_bed->report_info("Data read and FFT reports no error", 1);

  /** \todo Now test the resulting frequency is right.... Also tests our axes...*/
  
  //Get primary frequency
  
  test_bed->report_err(err);
  return err;

}

test_entity_basic_maths::test_entity_basic_maths(){

  name = "basic maths helpers";
  size = 256;
  data_square=(calc_type*)calloc(size,sizeof(calc_type));
  data_positive=(calc_type*)calloc(size,sizeof(calc_type));
  data_tmp=(calc_type*)calloc(size,sizeof(calc_type));
  axis=(calc_type*)calloc(size,sizeof(calc_type));
  d_axis=(calc_type*)calloc(size,sizeof(calc_type));

  axisf=(my_type*)calloc(size,sizeof(my_type));

  data_square[0] = 1.0;
  data_positive[0] = 0.0;
  axis[0] = 0.0;
  axisf[0] = 0.0;
  d_axis[0] = 1.0;
  for(int i=1; i<size; ++i){
    data_square[i] = - data_square[i-1];
    //alternating square wave, average=int=0
    data_positive[i] = data_positive[i-1] + 0.1;
    //monotonic increase, integral = size * size/10/2 depending on upper bnd
    d_axis[i] = 1.0;
    axis[i] = axis[i-1] + d_axis[i];
    axisf[i] = axisf[i-1] + 1.0;

  }

//set up some data arrays...

}
test_entity_basic_maths::~test_entity_basic_maths(){

  free(data_square);
  free(data_positive);
  free(data_tmp);
  free(axis);
  free(d_axis);
}

int test_entity_basic_maths::run(){
  int err=TEST_PASSED;

  my_type target;
  int whe;

  target = 13.5;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = 254.89;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = -2;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
  }
  target = 1302;
  whe = where(axisf, size, target);
  if(whe > 0){
    if(!(axisf[whe] >= target && axisf[whe-1] <= target)) err|=TEST_WRONG_RESULT;
    //test_bed->report_info(mk_str(target)+" "+mk_str(axisf[whe]), 2);
  }
  if(err == TEST_PASSED) test_bed->report_info("Where OK", 1);

  calc_type res = integrator(data_square, size, d_axis);
  if(res!= 0.0) err |= TEST_WRONG_RESULT;
  res = integrator(data_positive, size, d_axis);
  if(std::abs(res - (pow((calc_type)(size-1), 2)/20.0)) > res*PRECISION) err |= TEST_WRONG_RESULT;
  //test it's correct to within some finite precision, defined in header
  if(err == TEST_PASSED) test_bed->report_info("Integrator OK", 1);

  memcpy((void*)data_square, (void*)data_tmp, sizeof(calc_type)*size);

  inplace_boxcar_smooth(data_tmp, size, 2, 1);
  calc_type total=0;
  for(int i=0;i<size; ++i){
    total += data_tmp[i];
  }
  if(std::abs(total) > PRECISION) err |=TEST_WRONG_RESULT;
  //Smooth of 2 on square wave should give 0
  
  memcpy((void*)data_positive, (void*)data_tmp, sizeof(calc_type)*size);
  
  inplace_boxcar_smooth(data_tmp, size, 4, 0);
  total=0;
  for(int i=4;i<size-4; ++i){
    total += std::abs(data_positive[i] - data_tmp[i]);
  }
  //Smooth on straight line should do nothing except at ends...
  if(std::abs(total) > PRECISION) err |=TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Boxcar smooth OK", 1);


  {
    calc_type axis[4]={0.0,1.0,2.0,3.0}, vals[4]={0.0,1.0,0.0,1.0}, target, interp;

    target = 0.25;
    interp = interpolate(axis, vals, target, 1);
    if(std::abs(interp - 0.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 0.75;
    interp = interpolate(axis, vals, target, 1);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    //Nearest value interpol
    target = 0.5;
    interp = interpolate(axis, vals, target, 2);
    if(std::abs(interp - 0.5)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.0;
    interp = interpolate(axis, vals, target, 2);
    if(std::abs(interp - 1.0)> PRECISION) err |=TEST_WRONG_RESULT;
    target = 1.5;
    interp = interpolate(axis+1, vals+1, target, 2);
    if(std::abs(interp - 0.5)> PRECISION) err |=TEST_WRONG_RESULT;
    //2 pt linear interpol
  }
  
  if(err == TEST_PASSED) test_bed->report_info("Interpolate OK", 1);

  //  x^3 - 17x^2 + 92x - 150.
  std::vector<calc_type> result = cubic_solve(-17.0, 92.0, -150.0);
  if((result.size() != 1) || std::abs(result[0] - 3.0) > PRECISION) err|=TEST_WRONG_RESULT;
  //Test with polynomial of known integer coefficients and roots

  result = cubic_solve(-20.5, 100.0, -112.76);
  //test with random polynomial which happens to have 3 real roots
  calc_type res_el, tot;
  for(size_t i=0; i<result.size(); ++i){
    res_el = result[i];
    tot = std::pow(res_el, 3) -20.5*std::pow(res_el, 2) +100.0*res_el - 112.76;
    if(std::abs(tot) > PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Cubic root does not solve polynomial, mismatch "+mk_str(tot)+" for root "+mk_str(res_el), 2);
      
      test_bed->report_info("Vieta a) "+ mk_str(result[0]+result[1]+result[2] - 20.5), 2);
      test_bed->report_info("Vieta b) "+mk_str( result[0]*result[1]+result[0]*result[2]+result[2]*result[1] - 100.0), 2);
      test_bed->report_info("Vieta c) "+mk_str(result[0]*result[1]*result[2]-112.76), 2);

    }
  }
 	
  /*I think you should be able to recognize them using Vieta's formula for cubic equations, which states that if a cubic equation x3+ax2+bx+c=0x3+ax2+bx+c=0 has three different roots x1,x2,x3x1,x2,x3, then:
  −a=x1+x2+x3 b = x1x2+x1x3+x2x3 and -c = x1x2x3
**/
 
  test_bed->report_err(err);
  return err;

}

test_entity_extern_maths::test_entity_extern_maths(){

  name = "external maths routines";

}
test_entity_extern_maths::~test_entity_extern_maths(){

}

int test_entity_extern_maths::run(){

  //cyl_bessel_j(v, x) = Jv(x)
  //cyl_neumann(v, x) = Yv(x) = Nv(x)
  double bess, arg, bess1, bess2;
  int index, err=TEST_PASSED;

  index = 0;
  arg = 2.40482555769577;
  bess = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  index = 1;
  arg=7.01558666981561;
  bess = boost::math::cyl_bessel_j(index, arg);
  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  index=5;
  arg=12.3386041974669;
  bess = boost::math::cyl_bessel_j(index, arg);

  if(std::abs(bess - 0.0) >PRECISION) err|= TEST_WRONG_RESULT;

  //Test identities used to save time
  index=3;
  arg=1.3457174;
  bess = boost::math::cyl_bessel_j(index-1, arg);
  bess1 = boost::math::cyl_bessel_j(index+1, arg);
  bess2 = boost::math::cyl_bessel_j(index, arg);
  //test_bed->report_info(mk_str(std::abs(bess + bess1 - 2.0*(double)index*bess2/arg)), 1);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;

  index=7;
  arg=-20.98;
  bess = boost::math::cyl_bessel_j(index-1, arg);
  bess1 = boost::math::cyl_bessel_j(index+1, arg);
  bess2 = boost::math::cyl_bessel_j(index, arg);
  //test_bed->report_info(mk_str(std::abs(bess + bess1 - 2.0*(double)index*bess2/arg)), 1);
  if(std::abs(bess + bess1- 2.0*(double)index*bess2/arg) >PRECISION) err|= TEST_WRONG_RESULT;


  index =0;
  arg =1.0;
  bess = boost::math::cyl_bessel_j(index, arg);

  if(std::abs(bess - 0.7651976865579665514497) >PRECISION) err|= TEST_WRONG_RESULT;
  if(err == TEST_PASSED) test_bed->report_info("Bessel functions OK", 1);

  test_bed->report_err(err);

  return err;

}

test_entity_plasma::test_entity_plasma(){
  name = "plasma";
  plas = new plasma(-1.0);

}
test_entity_plasma::~test_entity_plasma(){

  delete plas;
}

int test_entity_plasma::run(){

  int err=TEST_PASSED;

  std::vector<calc_type> results;
  calc_type x=1.0, v_par, n=-1, om_ce_local, om_pe_local;
  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");
  
  v_par = 0.1* v0;
  calc_type cos_theta, mu_tmp1, mu_tmp2;
  cos_theta = std::cos(std::atan(x));
  calc_type gamma, gamma2;
  gamma2 = 1.0/( 1.0 - std::pow(v_par/v0, 2));
  
  gamma = std::sqrt(gamma2);
  
  results = plas->get_omega(x, v_par, n);
  /**Now check each element of this satisfies Stix 2.45 and the resonance condition together*/
  test_bed->report_info("Testing resonant frequency solver", 1);
  test_bed->report_info(mk_str((int)results.size())+" frequency solutions found", 2);
  
  for(int i=0; i<(int)results.size(); ++i){
    test_bed->report_info("Freq is "+mk_str(results[i])+" = "+mk_str(results[i]/my_const.omega_ce)+" om_ce", 2);
    
    mu_tmp1 = std::pow(v0 * (gamma*results[i] - n*om_ce_local)/(gamma*results[i] * v_par *cos_theta), 2);
    mu_tmp2 = (1.0 - (std::pow(om_pe_local,2)/(results[i]*(results[i] + om_ce_local*cos_theta))));
    if(std::abs(mu_tmp1 - mu_tmp2) > PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("refractive index mismatch of "+mk_str(mu_tmp1-mu_tmp2), 2);
    }
  }
  

  /** Now test if the returned mu matches the high density whistler in high dens regime */

  size_t n_tests = 10;
  calc_type tmp_omega=0.0, tmp_theta=pi/(calc_type)(n_tests), tmp_omega_n=0.0;
  mu_dmudom my_mu;
  mu my_mu_all;
  int err_cnt=0;
  test_bed->report_info("Testing whistler high density approx.", 1);

  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_ce_local)/(calc_type)(n_tests + 1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /** my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
      test_bed->report_info("Mismatch in high density approx or dispersion solver at "+mk_str(tmp_omega/std::abs(om_ce_local))+" "+mk_str(tmp_theta), 1);
      test_bed->report_info("Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu), 2);
    }
    //my_mu_all.mu and my_mu.mu should be exactly equal:
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }
  
  tmp_omega = 0.6*std::abs(om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_theta += pi/(calc_type)(n_tests);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /** my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
    
      test_bed->report_info("Mismatch in high density approx or dispersion solver at "+mk_str(tmp_omega/std::abs(om_ce_local))+" "+mk_str(tmp_theta), 1);
      test_bed->report_info("Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu), 2);
    }
     //my_mu_all.mu and my_mu.mu should be exactly equal:
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }

 
  }
  if(err_cnt> 0){
    test_bed->report_info("Total "+mk_str(err_cnt)+" out of "+mk_str(2*(int)n_tests)+" issues in high density approx or dispersion solver at precision: "+mk_str(LOW_PRECISION), 1);
    //err|=TEST_WRONG_RESULT;
    //Make these a warning not an error because we expect them sometimes
  }

  test_bed->report_info("Testing dispersion solver for plasma O mode", 1);

  //Try plasma wave modes in solvers, perpendicular propagation
  tmp_omega = om_pe_local;
  tmp_theta = pi/2.0;
  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_pe_local)/(calc_type)(n_tests + 1);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);

    mu_tmp2 = std::sqrt(std::pow(tmp_omega, 2) - std::pow(om_pe_local, 2))/tmp_omega;
    
    if(std::abs(my_mu_all.mu-mu_tmp2)/my_mu_all.mu > LOW_PRECISION){
      
      test_bed->report_info("Error in approx or dispersion solver for plasma wave at "+mk_str(tmp_omega/std::abs(om_pe_local))+" "+mk_str(tmp_theta), 1);
      
      test_bed->report_info("Mu "+mk_str(my_mu_all.mu)+" difference "+mk_str(my_mu_all.mu - mu_tmp2)+" relative error "+mk_str((my_mu_all.mu-mu_tmp2)/my_mu_all.mu), 2);
    }
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;

    }
    
  }
  //Try left hand X mode too
  test_bed->report_info("Testing dispersion solver for plasma X mode", 1);

  calc_type omega_UH = std::sqrt(om_pe_local*om_pe_local + om_ce_local*om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_pe_local)/(calc_type)(n_tests + 1);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta, false);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n, false);
    
    mu_tmp2 = std::sqrt(1.0 - std::pow(om_pe_local/tmp_omega, 2)*(std::pow(tmp_omega, 2) - std::pow(om_pe_local, 2))/(std::pow(tmp_omega, 2) - std::pow(omega_UH, 2)));
    if(std::abs(my_mu_all.mu-mu_tmp2)/my_mu_all.mu > LOW_PRECISION){
      
      test_bed->report_info("Error in approx or dispersion solver for plasma wave at "+mk_str(tmp_omega/std::abs(om_pe_local))+" "+mk_str(tmp_theta), 1);
      
      test_bed->report_info("Mu "+mk_str(my_mu_all.mu)+" difference "+mk_str(my_mu_all.mu - mu_tmp2)+" relative error "+mk_str((my_mu_all.mu-mu_tmp2)/my_mu_all.mu), 2);
    }
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;

    }
    
  }

/** \todo Test mu.dom and mu.dmudtheta too */
/** \todo Test phi? */



  test_bed->report_err(err);

  return err;

}

/**
For getting the resonant frequency to consider errors will matter, but the noise in spectra can be expected to be a similar sort of error
*
The reason for using the better dispersion solver is a) to avoid any numerical derivatives in general and b) because some of the calculation depends on derivs of the refractive index

*So basically I am using the easy approx where a) it’s a genuine nightmare to do better (10-14th order polynomial) and b) I expect other errors to be similar or more important and c) I really hope it’s linear or polynomial in error*/

test_entity_spectrum::test_entity_spectrum(){


}
test_entity_spectrum::~test_entity_spectrum(){


}

int test_entity_spectrum::run(){

  int err=TEST_PASSED;

  return err;

}

test_entity_albertG1::test_entity_albertG1(){

}
test_entity_albertG1::~test_entity_albertG1(){

}

int test_entity_albertG1::run(){
/** \todo WRITE!!!*/

  return TEST_PASSED;
}
