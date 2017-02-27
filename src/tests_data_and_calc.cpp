//
//  tests_data_and_calc.cpp
//  
//
//  Breakout by Heather Ratcliffe on 27/02/2017.
//
//

#ifdef RUN_TESTS_AND_EXIT

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <functional>
#include "tests_data_and_calc.h"
#include "reader.h"
#include "plasma.h"
#include "controller.h"
#include "data_array.h"
#include "spectrum.h"
#include "d_coeff.h"
#include "non_thermal.h"

#include <math.h>
#include <boost/math/special_functions.hpp>

test_entity_plasma::test_entity_plasma(){
  name = "plasma";
  plas = new plasma("./files/test");

}
test_entity_plasma::~test_entity_plasma(){

  delete plas;
}

int test_entity_plasma::run(){
/** \brief Test resonant frequencies and refractive indices
*
*Checks the resonant frequencies obey the equations used to derive them. Checks the dispersion roots for Whistlers match those found using high-density approx. Checks plasma O and X mode dispersion too. Note, first call for issues with these tests is to check returned mu.err on failing tests
*/

  int err=TEST_PASSED;
  
  err |= analytic_dispersion();
  err |= resonant_freq();
  err |= high_density();
  err |= other_modes();
  err |= phi_dom();

  return err;
}

int test_entity_plasma::analytic_dispersion(){
/** \brief Check analytic dispersion relations
*
*Checks the analytic relations, both ways and including derivatives \todo Check theta
*/

  int err=TEST_PASSED;

  calc_type k, om, om_new, d_om, d_k;
  size_t n_tests = 5;
  //First whistlers
  //At very large k we should get om_ce and at 0, we get 0
  k=100;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  if(std::abs(om /plas->get_omega_ref("ce") -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  k = 0;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  if(std::abs(om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  //Symmetry
  calc_type om_2;
  k=0.0001;
  om = plas->get_dispersion(k, WAVE_WHISTLER);
  k=-k;
  om_2 = plas->get_dispersion(k, WAVE_WHISTLER);
  
  if(std::abs(om - om_2) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  
  //If we pick a few e.g.s and do them both ways we should get same result back. Start with omega for simplicity
  for(size_t i=1; i<= n_tests; i++){
    om = (float)i *plas->get_omega_ref("ce")/(float) (n_tests +1);
    k = plas->get_dispersion(om, WAVE_WHISTLER, true);
    om_new = plas->get_dispersion(k, WAVE_WHISTLER);
    if(std::abs(om /om_new -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Derivs should be 0 at 0 and large omega.

  k=100;
  d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  k = 0;
  d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  
  //Finally derivs at our samples should be inverses
  for(size_t i=1; i<= n_tests; i++){
    om = (float)i *plas->get_omega_ref("ce")/(float) (n_tests +1);
    k = plas->get_dispersion(om, WAVE_WHISTLER, true);
    
    d_om = plas->get_dispersion(k, WAVE_WHISTLER, false, true);
    d_k = plas->get_dispersion(om, WAVE_WHISTLER, true, true);
    if(std::abs(d_om*d_k -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Finally we check how it handles out of range
  try{
    k = plas->get_dispersion(plas->get_omega_ref("ce")*1.5, WAVE_WHISTLER, true);
  }catch(...){
    err |= TEST_ASSERT_FAIL;
  }

//  std::cout<<"Err: "<<err<<'\n';
  //Now EM
  //At 0, we get om_pe
  k = 0;
  om = plas->get_dispersion(k, WAVE_O);
  if(std::abs(om - plas->get_omega_ref("pe")) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  
  //If we pick a few e.g.s and do them both ways we should get same result back. Go from say om_pe to 3om_pe
  for(size_t i=1; i<= n_tests; i++){
    om = plas->get_omega_ref("pe")*(1.0+ + 2.0*(float)i/(float) (n_tests +1));
    k = plas->get_dispersion(om, WAVE_O, true);
    om_new = plas->get_dispersion(k, WAVE_O);
    if(std::abs(om /om_new -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  
  //Derivs should be 0 at 0 and ~1/c at very large omega
  k = 0;
  d_om = plas->get_dispersion(k, WAVE_O, false, true);
  if(std::abs(d_om) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  om = 100* plas->get_omega_ref("pe");
  d_k = plas->get_dispersion(om, WAVE_O, true, true);
  if(std::abs(d_k*v0 - 1.0) > LOW_PRECISION) err |=TEST_WRONG_RESULT;
  std::cout<<"Err: "<<err<<'\n';

  //Finally derivs at our samples should be inverses
  for(size_t i=1; i<= n_tests; i++){
    om = plas->get_omega_ref("pe")*(1.0+ + 2.0*(float)i/(float) (n_tests +1));
    k = plas->get_dispersion(om, WAVE_O, true);
    
    d_om = plas->get_dispersion(k, WAVE_O, false, true);
    d_k = plas->get_dispersion(om, WAVE_O, true, true);
    if(std::abs(d_om*d_k -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  //Finally we check how it handles out of range
  try{
    k = plas->get_dispersion(plas->get_omega_ref("pe")*0.8, WAVE_O, true);
  }catch(...){
    err |= TEST_ASSERT_FAIL;
  }

  if(err == TEST_PASSED) test_bed->report_info("Analytic dispersion OK", 1);
  else  test_bed->report_info("Error in analytic dispersion", 1);

  return err;
}

int test_entity_plasma::resonant_freq(){
/** \brief Check resonant frequency solver
*
*Checks the returned resonant frequency obeys equations used to derive it by solving both for mu.
*/

  int err=TEST_PASSED;
  int n_tests = 5;
  std::vector<calc_type> results;
  mu_dmudom my_mu;
  int err_count = 0;
  calc_type x, v_par, n, om_ce_local, om_pe_local;
  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");


  calc_type cos_theta, mu_tmp1, mu_tmp2, tmp_omega_n=0;
  calc_type gamma, gamma2;

  test_bed->report_info("Testing resonant frequency solver", 1);


  for(int ii=0; ii<n_tests; ii++){
    v_par = (0.01 + 0.5*(float)ii/ (float)(n_tests+1))* v0;

    for(int j=0; j< n_tests; j++){
      x = 4.0* (float) j / (float)(n_tests+1);
      cos_theta = std::cos(std::atan(x));

      for(int k=0; k< n_tests; k++){
        n = -n_tests/2 + k*n_tests/2;
        
        gamma2 = 1.0/( 1.0 - std::pow(v_par/v0, 2));
        
        gamma = std::sqrt(gamma2);
        
        results = plas->get_resonant_omega(x, v_par, n);
        /**Now check each element of the resonant frequency solution set satisfies Stix 2.45 and the resonance condition together*/
        for(int i=0; i<(int)results.size(); ++i){
          //test_bed->report_info("Freq is "+mk_str(results[i], true)+" = "+mk_str(results[i]/my_const.omega_ce, true)+" om_ce", 2);
          
          mu_tmp1 = std::pow(v0 * (gamma*results[i] - n*om_ce_local)/(gamma*results[i] * v_par *cos_theta), 2);
          mu_tmp2 = (1.0 - (std::pow(om_pe_local,2)/(results[i]*(results[i] + om_ce_local*cos_theta))));

          if(std::abs((mu_tmp1 - mu_tmp2)/mu_tmp1) > NUM_PRECISION){
            err|=TEST_WRONG_RESULT;
            test_bed->report_info("refractive index mismatch of "+mk_str((mu_tmp1-mu_tmp2)/mu_tmp1), 2);
          }
        
          //Also check there is a valid full mu solution
          my_mu = plas->get_high_dens_phi_mu_om(results[i], std::atan(x), 0.0, 0.0, tmp_omega_n);
          if(my_mu.err){
            err|=TEST_WRONG_RESULT;
            test_bed->report_info("No full mu solution for resonant frequency", 2);
            err_count++;
          }
        }
      }
    }
  }
  test_bed->report_info("Mu error count: "+mk_str(err_count), 2);
  return err;

}

int test_entity_plasma::high_density(){
/** \brief Tests high density approximation for dispersion relations
*
*Test if the mu found by get_root and get_phi_mu_om matches the high density whistler in high dens regime
*/

  int err=TEST_PASSED;

  calc_type om_ce_local, om_pe_local;
  calc_type mu_tmp2;

  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");

  size_t n_tests = 10;
  calc_type tmp_omega=0.0, tmp_theta=pi/(calc_type)(n_tests+1), tmp_omega_n=0.0;
  mu_dmudom my_mu;
  mu my_mu_all;
  mu_dmudom my_mu_dens;
  int err_cnt=0;
  test_bed->report_info("Testing whistler high density approx.", 1);

  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_ce_local)/(calc_type)(n_tests + 1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_dens = plas->get_high_dens_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /** my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
      test_bed->report_info("Mismatch in high density approx or dispersion solver at omega="+mk_str(tmp_omega/std::abs(om_ce_local), true)+" om_ce and theta= "+mk_str(tmp_theta/pi)+" pi", 1);
      test_bed->report_info("Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu, true), 2);
    }
    /** my_mu_dens should EXACTLY equal Stix 2.45 without the 1.0 term*/
    mu_tmp2 = sqrt( - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));

    if(std::abs(my_mu_dens.mu-mu_tmp2)/my_mu_dens.mu > NUM_PRECISION){
      err_cnt++;
      test_bed->report_info("  Mismatch in high density approx or dispersion solver at omega="+mk_str(tmp_omega/std::abs(om_ce_local), true)+" om_ce and theta= "+mk_str(tmp_theta/pi)+" pi", 1);

      test_bed->report_info("  Mu "+mk_str(my_mu_dens.mu)+" difference "+mk_str(my_mu_dens.mu - mu_tmp2)+" relative error "+mk_str((my_mu_dens.mu-mu_tmp2)/my_mu_dens.mu), 2);
    }

    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.mu-my_mu.mu) > PRECISION){
      test_bed->report_info("Inconsistent root between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }
  
  tmp_omega = 0.6*std::abs(om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_theta += pi/(calc_type)(n_tests+1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);

    /* my_mu.mu should roughly equal Stix 2.45*/
    mu_tmp2 = sqrt(1.0 - (std::pow(om_pe_local,2)/(tmp_omega*(tmp_omega + om_ce_local*std::cos(tmp_theta)))));
    if(std::abs(my_mu.mu-mu_tmp2)/my_mu.mu > LOW_PRECISION){
      err_cnt++;
      test_bed->report_info("    Mismatch in high density approx or dispersion solver at omega="+mk_str(tmp_omega/std::abs(om_ce_local), true)+" om_ce and theta= "+mk_str(tmp_theta/pi)+" pi", 1);

      test_bed->report_info("    Mu "+mk_str(my_mu.mu)+" difference "+mk_str(my_mu.mu - mu_tmp2)+" relative error "+mk_str((my_mu.mu-mu_tmp2)/my_mu.mu), 2);
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

  return err;
}

int test_entity_plasma::other_modes(){
/** \brief Test dispersion solver with other wave modes */

  int err=TEST_PASSED;

  calc_type om_ce_local, om_pe_local;
  calc_type mu_tmp2;

  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");


  test_bed->report_info("Testing dispersion solver for plasma O mode", 1);
  size_t n_tests = 10;
  mu_dmudom my_mu;
  mu my_mu_all;

  /**Try plasma wave modes in solvers, perpendicular propagation*/
  calc_type tmp_omega = om_pe_local, tmp_omega_n=0;
  calc_type tmp_theta = pi/2.0;
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
  /**Try left hand X mode too*/
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

  return err;

}

int test_entity_plasma::phi_dom(){
/** \brief Test other plasma returns
*
* Checks the values of mu.dom, mu.dmudtheta and phi against special cases. \todo Write this Phi is pretty utestable really...
*/

  int err=TEST_PASSED;
  size_t n_tests = 10;

  calc_type mu_tmp1, mu_tmp2, om_ce_local, om_pe_local;
  om_ce_local = plas->get_omega_ref("ce");
  om_pe_local = plas->get_omega_ref("pe");

  calc_type d_omega = std::abs(om_ce_local)/1e8;
  calc_type d_theta = pi/(calc_type)(n_tests)/1e7;
  //Derivative step size.

  mu_dmudom my_mu, my_mu_p;
  mu my_mu_all, my_mu_all_p;

  calc_type tmp_omega = 0.0, tmp_theta=pi/(calc_type)(n_tests), tmp_omega_n=0.0;

  test_bed->report_info("Testing dmu/domega", 1);

  for(size_t i =0; i<n_tests; i++){
    tmp_omega += std::abs(om_ce_local)/(calc_type)(n_tests + 1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu_p = plas->get_phi_mu_om(tmp_omega+d_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all_p = plas->get_root(0.0, tmp_omega+d_omega, tmp_theta);

    /** Approx numerical derivative*/
    mu_tmp1 = (my_mu.mu - my_mu_p.mu)/d_omega;
    mu_tmp2 = (my_mu_all.mu - my_mu_all_p.mu)/d_omega;
    if(std::abs(std::abs(mu_tmp1/my_mu.dmudom) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_phi_mu_om", 2);
    }
    if(std::abs(std::abs(mu_tmp2/my_mu_all.dmudom) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_root", 2);
    }

    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.dmudom-my_mu.dmudom) > PRECISION){
      test_bed->report_info("Inconsistent derivative between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }
  
  test_bed->report_info("Testing dmu/dtheta", 1);
  tmp_theta = 0.0;
  tmp_omega = 0.7*std::abs(om_ce_local);
  for(size_t i =0; i<n_tests; i++){
    tmp_theta += pi/(calc_type)(n_tests+1);
    my_mu = plas->get_phi_mu_om(tmp_omega, tmp_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all = plas->get_root(0.0, tmp_omega, tmp_theta);
    my_mu_p = plas->get_phi_mu_om(tmp_omega, tmp_theta+d_theta, 0.0, 0.0, tmp_omega_n);
    my_mu_all_p = plas->get_root(0.0, tmp_omega, tmp_theta+d_theta);

    /** Approx numerical derivative. Manually fix signs*/
    mu_tmp1 = -(my_mu.mu - my_mu_p.mu)/d_theta;

    mu_tmp2 = -(my_mu_all.mu - my_mu_all_p.mu)/d_theta;
    
    if(std::abs(std::abs(mu_tmp1 /my_mu.dmudtheta) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_phi_mu_om at omega = "+mk_str(tmp_omega/std::abs(om_ce_local), true) +" and phi = "+mk_str(tmp_theta/pi, true)+" pi", 2);

    }
    if(std::abs(std::abs(mu_tmp2/my_mu_all.dmudtheta) - 1.0) > NUM_PRECISION){
      err|=TEST_WRONG_RESULT;
      test_bed->report_info("Wrong derivative in get_root at omega = "+mk_str(tmp_omega/std::abs(om_ce_local), true) +" and phi = "+mk_str(tmp_theta/pi, true)+" pi", 2);

    }
    /**my_mu_all.mu and my_mu.mu should be exactly equal*/
    if(std::abs(my_mu_all.dmudtheta-my_mu.dmudtheta) > PRECISION){
      test_bed->report_info("Inconsistent derivative between get_root and get_phi_mu_om", 2);
      err|=TEST_WRONG_RESULT;
    }
    
  }

  return err;

}

/*
For getting the resonant frequency to consider errors will matter, but the noise in spectra can be expected to be a similar sort of error
*
The reason for using the better dispersion solver is a) to avoid any numerical derivatives in general and b) because some of the calculation depends on derivs of the refractive index

*So basically I am using the easy approx where a) it’s a genuine nightmare to do better (10-14th order polynomial) and b) I expect other errors to be similar or more important and c) I really hope it’s linear or polynomial in error*/
//----------------------------------------------------------------

test_entity_spectrum::test_entity_spectrum(){

  name = "spectrum checks";
  file_prefix = "./files/";
  
}
test_entity_spectrum::~test_entity_spectrum(){

  if(test_contr) delete test_contr;

}

int test_entity_spectrum::run(){
/** \brief Test spectrum extraction
*
* This should test the dispersion relation approximations are OK (plain and vg). Check test spectrum makes sense. Test extraction of a spectrum from data. Note data does not come from files, but from a test file which is already written as a data array using ./files/generate_fftd.pro which makes FFT_data.dat with the FFTd data and spectrum.dat with a derived spectrum to check against.
*/

  int err = TEST_PASSED;

  //Use a different deck.status file...
  deck_constants const_tmp = my_const;
  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
  share_consts();

  err|= setup();
  if(!test_bed->check_for_abort(err)) err|= basic_tests1();
  if(!test_bed->check_for_abort(err)) err|= basic_tests2();
  if(!test_bed->check_for_abort(err)) err|= albertGs_tests();
  
  my_const = const_tmp;
  share_consts();

  return err;

}

int test_entity_spectrum::setup(){
/** \brief Setup to test spectrum
*
* Note strictly this is the test of data array constructor taking a filename too.
*/

  int err = TEST_PASSED;

  test_dat_fft = data_array(file_prefix + "FFT_data.dat", true);
  if(!test_dat_fft.is_good()){
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }
  test_contr = new controller(file_prefix);

  return err;
}

int test_entity_spectrum::basic_tests1(){
/** \brief Basic tests of spectrum
*
*Test test_spectrum and angle generation \todo The angles are integrating to 0.5 not 1. Which do we want????
*/
  int err = TEST_PASSED;

  std::fstream outfile;
  size_t len=0;
  my_type total_error =0.0;
  my_type * d_angle, * angle_data;

  test_contr->add_spectrum(test_dat_fft.get_dims(0), DEFAULT_N_ANG, true);
  //Delete and re-add to test that
  test_contr->delete_current_spectrum();
  test_contr->add_spectrum(test_dat_fft.get_dims(0), DEFAULT_N_ANG, true);

  /** Check this test spectrum makes sense \todo HOW????*/

  test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_DELTA);
  //Check angle distrib integrates to 1 for each case
  //NOTE we can only do this if MIN_ANG is either 0 or is - MAX_ANG. otherwise we're into erf and bunk
  bool is_symmetric=false, is_zero = false;
  if(std::abs(ANG_MIN + ANG_MAX) < PRECISION) is_symmetric = true;
  if(std::abs(ANG_MIN) < PRECISION) is_zero = true;
  
  if(is_symmetric || is_zero){
    len = test_contr->get_current_spectrum()->get_g_dims(1);
    d_angle = (my_type *) calloc(DEFAULT_N_ANG, sizeof(my_type));
    for(size_t i=0; i<DEFAULT_N_ANG-1; ++i){
      d_angle[i] = std::abs(test_contr->get_current_spectrum()->get_ang_axis_element(i) - test_contr->get_current_spectrum()->get_ang_axis_element(i+1));
    }
    angle_data = (my_type *) malloc(len*sizeof(my_type));
    for(size_t i=0; i<len; i++){
      *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);
    }
    
    total_error = integrator(angle_data, len, d_angle);
    test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_GAUSS);
    for(size_t i=0; i<len; i++) *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);

    total_error += integrator(angle_data, len, d_angle);

    test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_ISO);
    for(size_t i=0; i<len; i++) *(angle_data + i) = test_contr->get_current_spectrum()->get_g_element(0, i);
    total_error += integrator(angle_data, len, d_angle);
    
    my_type expected = is_zero ? 2.0 : 3.0;
    //Iso always integrates to 1. Gaussian and delta are always symmetric
    if(std::abs(total_error - expected)/3.0 > NUM_PRECISION){
    
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("Error in angular distribution integrals, value " + mk_str(total_error, true));
    }
    if(angle_data) free(angle_data);
    if(d_angle) free(d_angle);
  }else{
    test_bed->report_info("Cannot test assymmetric spectrum");
  }
  outfile.open("spect_testy.dat", std::ios::binary|std::ios::out);
  test_contr->get_current_spectrum()->write_to_file(outfile);
  outfile.close();
  if(err == TEST_PASSED) test_bed->report_info("Test spectrum OK");
  
  /** Now make the real spectrum from data and check the result matches the plain text test file*/
  return err;
}
int test_entity_spectrum::basic_tests2(){
/**\brief Test spectrum extraction
*
*Compare extracted spectrum from FFT'd data file to test file
*/
  int err = TEST_PASSED;

  std::fstream outfile, infile;
  size_t len=0;
  my_type total_error =0.0;

  test_contr->get_current_spectrum()->generate_spectrum(test_dat_fft ,10, FUNCTION_GAUSS);

  test_spect = data_array(file_prefix + "spectrum.dat", true);
  if(test_spect.is_good()){
    //We ignore frequencies below say 0.05 om_ce
    my_type * ax = test_spect.get_axis(0, len);
    int min_ind = 0;
    if(ax) min_ind = where(ax+len/2, len/2, 17588.200*0.05);
    /**Hard code min freq to match the IDL file with test data generation...*/
    
    total_error = 0.0;
    for(size_t i=0; i< len/2 - min_ind; i++){
      total_error += std::abs(test_contr->get_current_spectrum()->get_B_element(i)-test_spect.get_element(i));
    }
    for(size_t i=len/2 + min_ind; i< len; i++){
      total_error += std::abs(test_contr->get_current_spectrum()->get_B_element(i)-test_spect.get_element(i));

    }
    if(total_error > LOW_PRECISION){
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("Mismatch between generated spectrum and test spectrum of "+mk_str(total_error));
    }
    /* Preserve the spectrum*/
    outfile.open("spect_out.dat", std::ios::binary|std::ios::out);
    test_contr->get_current_spectrum()->write_to_file(outfile);
    outfile.close();
  }else{
    err |= TEST_ASSERT_FAIL;
  }
  if(err == TEST_PASSED) test_bed->report_info("Generate spectrum OK");
  
  data_array old_B = test_contr->get_current_spectrum()->copy_out_B();
  //Now dump to file and read back in and compare
  bool err2 = test_contr->add_spectrum("spect_out.dat");
  if(err2) err |= TEST_ASSERT_FAIL;
  else{
    data_array new_B = test_contr->get_current_spectrum()->copy_out_B();
    if(!old_B.is_good() || !new_B.is_good() || (old_B != new_B)){
      test_bed->report_info("Error or Mismatch in read", 0);
      err|= TEST_WRONG_RESULT;
    }
  }
  return err;

}
int test_entity_spectrum::technical_tests(){

}
int test_entity_spectrum::albertGs_tests(){
/** \brief Tests of the Albert G functions in spectrum. Also tests the normalisations on the way. NOTE: since we're comparing the values of an analytic function with a numerical integral, we can get mismatches at the cutoffs. More points should help this. If that doesn't there may be something wrong.
*
*
*/
  int err = TEST_PASSED;

  calc_type om_ce_local, om_pe_local, G1, G2, G1_analytic, G2_analytic;
  om_ce_local = test_contr->get_plasma().get_omega_ref("ce");
  om_pe_local = test_contr->get_plasma().get_omega_ref("pe");

  calc_type mass_ratio = 1836.2;

  size_t n_tests = 10;
  calc_type tmp_omega=0.0, tmp_x;

  test_contr->add_spectrum(5000, DEFAULT_N_ANG, true);

  if(!test_contr->get_current_spectrum()->is_good()){
    my_error_print("Spectrum in invalid state. Aborting", mpi_info.rank);
    err |=TEST_ASSERT_FAIL;
    err |=TEST_FATAL_ERR;
    return err;
  }

  test_contr->get_current_spectrum()->make_test_spectrum(FUNCTION_GAUSS);
  
  my_type om_min, om_max, x_min, x_max, om_peak;
  om_min = 12000.0;
  om_max = 16500.0;
  //make sure this is lower than the test spectrum axis range
  x_min = 0.0;
  x_max = 3.0;
  
  test_contr->get_current_spectrum()->truncate_om(om_min, om_max);
  test_contr->get_current_spectrum()->truncate_x(x_min, x_max);
  om_peak = test_contr->get_current_spectrum()->get_peak_omega();
  //Now we have a test spectrum. Need to know what its normalisations should be. And what the Albert functions should resolve to.
  
  my_type width=0.1*om_peak;
  
  for(size_t i=0; i< n_tests;i++){
    tmp_omega = std::abs(om_ce_local)/10.0 + 89.0/100.0 * om_max * (1.0 - exp(-i));
    //Cover range from small to just below om_ce...
    G1 = test_contr->get_current_spectrum()->get_G1(tmp_omega);

    //Analytic calculations for truncated Gaussians, see Albert 2005
    
    if(tmp_omega > om_min && tmp_omega < om_max){
      G1_analytic = 2.0 / std::sqrt(pi) * std::exp( - std::pow((tmp_omega - om_peak)/width, 2));
      G1_analytic /= (boost::math::erf((om_max - om_peak)/width) +boost::math::erf((om_peak - om_min)/width));
      G1_analytic /=width;
    }else{
      G1_analytic = 0.0;
    }

    if( (G1 != 0.0 && std::abs(G1-G1_analytic)/(G1) > LOW_PRECISION)|| (G1 == 0.0 && G1_analytic != 0.0)){
      err |= TEST_WRONG_RESULT;
      test_bed->report_info("G1 does not match analytic calc, relative error = "+mk_str((std::abs(G1/G1_analytic)-1.0)*100, true)+" at "+mk_str(tmp_omega, true), mpi_info.rank);
    }
  }


  tmp_omega = 0.6 * std::abs(om_ce_local);
  for(size_t i=0; i< n_tests;i++){

    tmp_x = ANG_MIN + i * (ANG_MAX - ANG_MIN)/(n_tests-1);
    
    G2 = test_contr->get_current_spectrum()->get_G2(tmp_omega, tmp_x);
    
  //  std::cout<<"  G2 "<<G2<<std::endl;
    G2_analytic = std::pow((( mass_ratio / (1.0 + mass_ratio))*om_ce_local*om_ce_local/om_pe_local/om_pe_local), 1.5);
    G2_analytic *= exp(- (tmp_x*tmp_x)/std::pow(SPECTRUM_ANG_STDDEV, 2));
    
    }


  std::fstream outfile;
  outfile.open("spect_truncated.dat", std::ios::out|std::ios::binary);
  test_contr->get_current_spectrum()->write_to_file(outfile);
  outfile.close();


  return err;
}
//----------------------------------------------------------------

test_entity_levelone::test_entity_levelone(){

  name = "level-one derivation";
  
  test_contr=nullptr;
  my_reader=nullptr;
  
}
test_entity_levelone::~test_entity_levelone(){

  delete test_contr;
  delete my_reader;
}

int test_entity_levelone::run(){
/** \brief Test entire level-1 data extraction
*
**/

  int err = TEST_PASSED;

  //Use a different deck.status file...
//  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
//  share_consts();

  strcpy(block_id, "ay");

  file_prefix = "./files/l1/l1";
  space_in[0] = 0;
  space_in[1] = 1024;
  time_in[0] = 0;
  time_in[1] = 4;
  time_in[2] = 100;

  err|= setup();
  if(!test_bed->check_for_abort(err)) err|= basic_tests();
  if(my_reader) delete my_reader;
  if(test_contr) delete test_contr;
    
  if(!test_bed->check_for_abort(err)){
    file_prefix = "./files/2dtest/";
    strcpy(block_id, "ey");
    space_in[0] = 0;
    space_in[1] = 1024;
    time_in[0] = 0;
    time_in[1] = 50;
    time_in[2] = 0;

    err|=setup();
    if(!test_bed->check_for_abort(err)) err|= twod_tests();
    if(!test_bed->check_for_abort(err)) err|= twod_space_tests();

  }
  
  return err;
}

int test_entity_levelone::setup(){
/** \brief Setup to "level one" extraction
*
*
*/

  int err = TEST_PASSED;
  bool use_row_time=false;
  my_reader = new reader(file_prefix, block_id);
  if(my_reader->current_block_is_accum()) use_row_time = true;

  if(!use_row_time){
    n_tims = std::max((int) (time_in[1]-time_in[0]), 1);
  }else{
    n_tims = time_in[2];
  }

  size_t n_dims;
  std::vector<size_t> dims;
  int err2 = my_reader->read_dims(n_dims, dims);
  if(err2) err |= TEST_FATAL_ERR;
  
  test_contr = new controller(file_prefix);

  return err;
}

int test_entity_levelone::basic_tests(){
/** \brief Basic tests of process to make levl-1 data
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  int space_dim = space_in[1]-space_in[0];

  dat = data_array(space_dim, n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_error_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

  int err2 = my_reader->read_data(dat, time_in, space_in);
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array...
  dat_fft = data_array(space_dim, n_tims);

  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }
  dat.B_ref = get_ref_Bx(file_prefix, space_in, 1);
  //Use first file. 0th accumulator is perhaps broken
  err2 = fft_array(dat,dat_fft);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->set_plasma_B0(dat.B_ref);
  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->generate_spectrum(dat_fft);
  
  int n_dims = dat.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=3){
    lims.push_back(-0.002);
    lims.push_back(0.002);
  }
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);
  
  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
//    dat_fft.write_section_to_file(file, lims);
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}

int test_entity_levelone::twod_tests(){
/** \brief Basic tests of process to make levl-1 data from 2-d input
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  size_t n_dims_in;
  std::vector<size_t> dims_in;
  my_reader->read_dims(n_dims_in, dims_in);
  if(n_dims_in != 2){
    test_bed->report_info("Wrong file dimension", 1);
    return TEST_FATAL_ERR;
  }
  int space_dim = space_in[1]-space_in[0];
//  dat = data_array(space_dim, dims_in[1], n_tims);
  dat = data_array(space_dim, n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_error_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

//  int err2 = my_reader->read_data(dat, time_in, space_in);
  int err2 = my_reader->read_data(dat, time_in, space_in, 1);
  //Note this totals on y-dim automagically
  
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array and size FFT accordingly
//  dat_fft = data_array(space_dim, dims_in[1], n_tims);
  dat_fft = data_array();
  dat_fft.clone_empty(dat);
  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }
  dat.B_ref = -1;
  err2 = fft_array(dat,dat_fft);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->make_test_spectrum();
  
  int n_dims = dat.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=3){
    lims.push_back(-0.002);
    lims.push_back(0.002);
  }
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);

  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}
int test_entity_levelone::twod_space_tests(){
/** \brief Basic tests of process to make levl-1 data from 2-d input
*
* Reads proper data files, produces FFT, derived spectrum etc*/
  int err = TEST_PASSED;

  size_t n_dims_in;
  std::vector<size_t> dims_in;
  my_reader->read_dims(n_dims_in, dims_in);
  if(n_dims_in != 2){
    test_bed->report_info("Wrong file dimension", 1);
    return TEST_FATAL_ERR;
  }
  int space_dim = space_in[1]-space_in[0];
  dat = data_array(space_dim, dims_in[1], n_tims);
  strcpy(dat.block_id, block_id);

  if(!dat.is_good()){
    my_error_print("Data array allocation failed.", mpi_info.rank);
    err |= TEST_ASSERT_FAIL;
    err |= TEST_FATAL_ERR;
  }

  int err2 = my_reader->read_data(dat, time_in, space_in);
  if(err2 == 1){
    return TEST_FATAL_ERR;
  }
  if(err2 == 2) n_tims = dat.get_dims(1);
  //Check if we had to truncate data array and size FFT accordingly
//  dat_fft = data_array(space_dim, dims_in[1], n_tims);

  dat.B_ref = -1;

//  dat = dat.total(2);

  dat_fft.clone_empty(dat);
  if(!dat_fft.is_good()){
    return TEST_FATAL_ERR;
  }

  err2 = fft_array(dat,dat_fft);

  dat_fft = dat_fft.total(2, 0.01f*my_const.omega_ce, 1.5f*my_const.omega_ce);

  test_bed->report_info("FFT returned err_state " + mk_str(err2));

  test_contr->add_spectrum(space_dim, DEFAULT_N_ANG, true);
  test_contr->get_current_spectrum()->make_test_spectrum();
  
  size_t n_dims = dat_fft.get_dims();
  std::vector<my_type> lims;
  if(n_dims >=2){
    lims.push_back(-0.002);
    lims.push_back(0.002);
    lims.push_back(-0.002);
    lims.push_back(0.002);

  }
  if(n_dims >=3){
    lims.push_back(-3.0*my_const.omega_ce);
    lims.push_back(3.0*my_const.omega_ce);
  }
  
//Set cutout limits on FFT
  std::string filename, time_str;
  time_str = mk_str(dat_fft.time[0], true)+"_"+mk_str(this->n_tims);
  std::string block = block_id;
  filename = file_prefix+"FFT_k_"+block +"_"+time_str+"_"+mk_str(dat_fft.space[0])+"_"+mk_str(dat_fft.space[1]) + ".dat";
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::binary);
  if(file.is_open()){
    dat_fft.write_section_to_file(file, lims);
//    dat.write_to_file(file);
    if(err2){
      test_bed->report_info("File writing failed");
      err |=TEST_ASSERT_FAIL;
    }
    
  }else{
    err |=TEST_ASSERT_FAIL;
  
  }
  file.close();
  test_bed->report_info("FFT section output in "+filename, 1);

  return err;

}

//----------------------------------------------------------------

test_entity_d::test_entity_d(){

  name = "D checks";
  file_prefix = "./files/";

}
test_entity_d::~test_entity_d(){

  if(test_contr) delete test_contr;
}

int test_entity_d::run(){
/** \todo WRITE!*/
  int err = TEST_PASSED;
  
  deck_constants const_tmp = my_const;
  if(mpi_info.rank == 0) get_deck_constants(file_prefix);
  share_consts();

  test_contr = new controller(file_prefix);
  test_bed->report_info("Reading spectrum", mpi_info.rank);
  //Now dump to file and read back in and compare
  bool err2 = test_contr->add_spectrum("spect_out.dat");
  if(!err2){
    test_bed->report_info("Calculating test D", mpi_info.rank);
    test_contr->add_d(5, 5);
    d_report report = test_contr->get_current_d()->calculate();
    err |= report.error;
    test_bed->report_info("Writing test file", mpi_info.rank);

    std::fstream file;
    file.open("test_d.dat", std::ios::binary|std::ios::trunc|std::ios::out|std::ios::in);
    if(file) test_contr->get_current_d()->write_to_file(file);
    file.close();
  }else{
    err |= TEST_ASSERT_FAIL;
  }
  my_const = const_tmp;
  share_consts();

  return err;
}

//----------------------------------------------------------------

test_entity_bounce::test_entity_bounce(){

}
test_entity_bounce::~test_entity_bounce(){

}

int test_entity_bounce::run(){
/** \todo write!*/
  int err = TEST_PASSED;

  return err;

}

test_entity_nonthermal::test_entity_nonthermal(){
  name = "nonthermal";
}

test_entity_nonthermal::~test_entity_nonthermal(){

}

int test_entity_nonthermal::run(){
/** \brief Check nonthermal module
*
* Check function binding, df/dp values, lookup table
*/
  int err = TEST_PASSED;
  
  //Setup single Max non-thermal
  non_thermal * my_elec = new non_thermal("./files/test");

  calc_type a_par = std::sqrt(2.0)*my_elec->get_v_par() / std::sqrt(1.0 - std::pow(my_elec->get_v_par()/v0, 2));
  calc_type a_perp_sq = 2.0*std::pow(my_elec->get_v_perp(), 2)/(1.0 - 2.0*std::pow(my_elec->get_v_perp()/v0, 2));

  calc_type norm = (my_elec->get_total_dens() - 1.0)/pi/std::sqrt(pi)/a_par/a_perp_sq;

  calc_type df_tmp, f_tmp = my_elec->f_p(0, 0);
  
  if(std::abs(f_tmp / norm -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  const size_t n_trials = 5;
  calc_type inc = 0.1, p_par, p_perp;

  //Test f against expected values
  for(size_t i=0; i<n_trials*n_trials; i++){
    p_perp = inc* (i % n_trials);
    p_par = inc* (int)(i/n_trials);
    f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
    if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }
  //Test a couple of extreme cases
  p_perp = 1.5;
  f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
  if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  p_par = 2.5;
  f_tmp = my_elec->f_p(p_par*a_par, p_perp*std::sqrt(a_perp_sq));
  if(std::abs(f_tmp /( norm*std::exp(-p_par*p_par)*std::exp(-p_perp*p_perp)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

  //Test f against expected values
  for(size_t i=0; i<n_trials; i++){
    p_perp = inc*i + inc;
    p_par = inc*i + inc;
    //close to 0 deriv -> 0

    f_tmp = my_elec->d_f_p(p_par*a_par, 0, 1);
    if(std::abs(f_tmp /((-2.0)*p_par/a_par*my_elec->f_p(p_par*a_par, 0)) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;

    f_tmp = my_elec->d_f_p(0, p_perp*std::sqrt(a_perp_sq), 0);
    if(std::abs(f_tmp /((-2.0)*p_perp/std::sqrt(a_perp_sq)*my_elec->f_p(0, p_perp*std::sqrt(a_perp_sq))) -1.0) > GEN_PRECISION) err |=TEST_WRONG_RESULT;
  }

  //This next is deep wings so we relax constraint so we're not trying to resolve 1e-60
  df_tmp = my_elec->d_f_p(0, 10.0*std::sqrt(a_perp_sq), 0);
  if(std::abs(df_tmp /((-2.0)*10.0/std::sqrt(a_perp_sq)*my_elec->f_p(0, 10.0*std::sqrt(a_perp_sq))) -1.0) > 2e-6) err |=TEST_WRONG_RESULT;
  
  if(err == TEST_PASSED) test_bed->report_info("Nonthermal values OK",1);
  
  err |= test_lookup();
  return err;
}

int test_entity_nonthermal::test_lookup(){
  /** Check known lookup function against known analytic results. Data generated by generate_lookup in generate_test_data.pro*/
  int err = TEST_PASSED;
  
  non_thermal * my_elec_an = new non_thermal("./files/an");
  non_thermal * my_elec_lookup = new non_thermal("./files/lookup");

  calc_type f_an, f_lookup, p_par, p_perp;
  const size_t n_tests = 100;

  for(size_t i=0; i< n_tests; i++){
    for(size_t j=0; j< n_tests; j++){
      p_par = (float) i *v0/2.0/(float) n_tests;
      p_perp = (float) j *v0/2.0/(float) n_tests;
      f_an = my_elec_an->f_p(p_par, p_perp);
      f_lookup = my_elec_lookup->f_p(p_par, p_perp);
      if(std::abs(f_lookup/f_an - 1.0) > LOW_PRECISION){
        err |= TEST_WRONG_RESULT;
      }
    }
  }
  if(err == TEST_PASSED) test_bed->report_info("Lookup OK", 1);
  return err;
}

#endif