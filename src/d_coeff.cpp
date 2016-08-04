//
//  d_coeff.cpp
//  
//

#include <math.h>
#include <algorithm>
#include <boost/math/special_functions.hpp>
#include "support.h"
#include "d_coeff.h"

extern const mpi_info_struc mpi_info;

extern deck_constants my_const;

diffusion_coeff::diffusion_coeff(int nx, int n_angs):data_array(nx, n_angs){
/** \brief Create coefficient
*
*Creates rectangular data_array and sets additional parameters
*/
  my_controller = nullptr;

  n_thetas = 100;
  n_n = 10;
  tag = LOCAL;
  // FAKENUMBERS for testing
}

diffusion_coeff::~diffusion_coeff(){


}

void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE]){
/**\brief Set parameters
*
*Sets the time and space ranges, wave type etc attached to the spectrum. Times should be in terms of file output time. Space in terms of grid points.
*/

  this->time[0] = time1;
  this->time[1] = time2;
  this->space[0] = space1;
  this->space[1] = space2;
  strcpy(this->block_id, block_id);
  this->wave_id = wave_id;
}


bool diffusion_coeff::write_to_file(std::fstream &file){
/** \todo Add any additionals */

  if(!file.is_open()) return 1;
  data_array::write_to_file(file);
  return 0;

}

void diffusion_coeff::make_velocity_axis(){
/**\brief Set velocity axis
*
*Makes linear axis between V_MIN and V_MAX
*/
  calc_type res = (V_MAX - V_MIN)/this->get_dims(0);
  long offset = std::abs(V_MIN)/res;
  make_linear_axis(0, res, offset);
}

void diffusion_coeff::make_pitch_axis(){
/**\brief Set pitch angle axis (tan theta)
*
*Makes linear axis between ANG_MIN and ANG_MAX
*/

  calc_type res = (ANG_MAX - ANG_MIN)/this->get_dims(1); //To cover range from Min to Max
  long offset = ANG_MIN/res;
  make_linear_axis(1, res, offset);
}

void diffusion_coeff::calculate(){
/** \brief Calculate D from wave spectrum and plasma
*
*Uses the data available via my_controller to calculate D, the raw diffusion coefficient as function of particle velocity and pitch angle.
* \todo Can we break this down at all? \todo n ranges
*/

/*Subsections:
To solve Eq 2 in Albert:
Get G1

Get phi (note also needs mu internally)
Get mu, dmu/domega which are used to:
  get G2 - needs mu, dmu/domega
  get denominator | 1 - (d omega / d kparallel)_x | from 6 - needs mu, dmu/domega
  get Eq 7 factor - needs mu, dmu/domega, dmu/dx

*/

//loop over th, w, psi, alpha, n
//HERE we want to end up with D_alpha_alpha

  //To DO!
  //get n ranges... and add in

  //We probably want to set the n and theta ranges in another routine... or using UI data

  plasma * plas;
  spectrum * spect;
  if(my_controller){
    plas = my_controller->get_plasma();
    spect = my_controller->get_current_spectrum();
  }
  else{
    my_print("No controller", mpi_info.rank);
    return;
  }

  this->copy_ids(spect);
  //copy block id, ranges etc from spect.
  calc_type theta, lat, omega_n=0.0, D_tmp, k_thresh;
  calc_type alpha, v_par, c2th, s2alpha; /** temporaries for clarity*/
  calc_type Eq6, mu_dom_mu, Eq7, dmudx, numerator;

  int n_min, n_max;
  std::vector<calc_type> omega_calc;
  mu_dmudom my_mu;
  
  calc_type om_ce_ref = plas->get_omega_ref("ce");
  calc_type D_consts = 0.5* pi*om_ce_ref*pow(v0, 3);
  //Time saving constant

  //Temporaries for wave normal angle***********************************
  calc_type *D_theta = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  calc_type *x = (calc_type *) calloc(n_thetas, sizeof(calc_type));
 
  /** Why are we using different angles for D and for waves?*/

  for(int i=0; i<n_thetas; ++i){
    x[i] = i* 4.0/n_thetas;
  }

  calc_type *dx = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  for(int i=0; i<n_thetas -1; ++i){
    dx[i] = x[i+1] - x[i];
  }
  
  
  k_thresh = spect->check_upper();
  //*******************************************************************

  
  //innermost loop should be n. Next theta, as we need mu at each theta. Omega and x are interchangeable...
  //Alpha remains, as does particle v.


  int last_report=0;
  int report_interval = dims[0]/10; //10 prints per round
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;

  int counter = 0, non_counter = 0;

//-------------------Main loops here----------------------------
//We have deep nested loops. Move ANYTHING that can be as far up tree as possible

//  for(int i =0; i< ((1< dims[0]) ? 1:dims[0]); ++i){
  for(size_t i =0; i< dims[0]; ++i){
    //particle parallel velocity
    v_par = get_axis_element(0, i);
    //Get limits on n for each velocity
    n_min = get_min_n(v_par, k_thresh, om_ce_ref);
    n_max = get_max_n(v_par, k_thresh, om_ce_ref);
    my_print("Velocity "+mk_str(v_par/v0, true)+" c", mpi_info.rank);

    if((i-last_report) >= report_interval){
      my_print("i "+mk_str(i), mpi_info.rank);
  
      last_report = i;
    }

//    for(int k =0; k< ((1 <dims[1]) ? 1: dims[1]); k++){
    for(size_t k =0; k< dims[1]; k++){
      //particle pitch angle
      alpha = get_axis_element(1, k);
      s2alpha = std::pow(std::sin(alpha), 2);

      for(int j=0;j<n_thetas; ++j){
      //theta loop for wave angle or x=tan theta
        theta = atan(x[j]);
        c2th = std::pow(cos(theta), 2);

        D_tmp = 0.0;
        for(int n=n_min; n<n_max; ++n){
          // n is resonant number
          omega_calc = plas->get_resonant_omega(x[j], v_par, (calc_type) n);

          for(size_t ii =0; ii< omega_calc.size(); ++ii){
          //each solution
            //std::cout<<"Freq is "<<omega_calc[ii]/my_const.omega_ce<<std::endl;

            my_mu = plas->get_high_dens_phi_mu_om(omega_calc[ii], theta, alpha, n, omega_n);
            if(my_mu.err){
              //Once angle is included we have no solution
              non_counter++;
              continue;
              
            }else{
              counter ++;
            }
            //No solution
            
            mu_dom_mu = my_mu.mu + omega_calc[ii] * my_mu.dmudom;
            dmudx = my_mu.dmudtheta *c2th;
            //Chain rule...

            // FAKENUMBERS
            Eq6 = omega_calc[ii]/(omega_calc[ii] - omega_n)* my_mu.mu/mu_dom_mu;

            Eq7 = -1.0* (my_mu.mu*omega_n/(omega_calc[ii]*(omega_calc[ii]-omega_n)) - my_mu.dmudom)/(my_mu.mu *std::sin(theta)*std::cos(theta) - dmudx);
            //Need this iff we use second expression in Eq 5
            numerator = std::pow( -s2alpha + omega_n/omega_calc[ii], 2);
            D_tmp += numerator * my_mu.phi/std::abs(1.0 - Eq6)*spect->get_G1(omega_calc[ii])*spect->get_G2(omega_calc[ii], x[j]); // FAKENUMBERS This will be the n summed D so add onto it each n iteration
          }
        }
        //Store into theta array. Might not need tmp in the end
        //Except that it saves us one multiplication per iteration without acting on element in place after.
        D_theta[j] = D_tmp*c2th;

      }
      //now integrate in x = tan theta
      D_tmp = integrator(D_theta, dims[1], dx);
      //Leaves us with one D for each alpha and v.
      //Store into data at i, k
      set_element(i, k, D_tmp*D_consts);
      //What about 1/delta omega???
    
    }
  }
//------------End main loops-----------------------------------

  free(dx);
  free(x);
  free(D_theta);
  
  my_print(mk_str(counter) + " solutions vs "+ mk_str(non_counter), mpi_info.rank);
  
}

int diffusion_coeff::get_min_n(calc_type v_par, my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs) \todo Is there a tighter limt? This is quite weak...
*/

  calc_type gamma = gamma_rel(v_par);
  /** \todo FIX!!! */

  return std::max(-(int)(gamma * k_thresh * std::abs(v_par / om_ce)), -n_n);

}

int diffusion_coeff::get_max_n(calc_type v_par, my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs)
*/

  calc_type gamma = gamma_rel(v_par);
  /** \todo FIX!!! */

  return std::min((int)(gamma * k_thresh * std::abs(v_par / om_ce)), n_n);


}

void diffusion_coeff::copy_ids( spectrum * src){
/** Copies ID fields from src array to this*/

  strcpy(this->block_id, src->block_id);
  std::copy(src->time, src->time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src->space[i];
}


